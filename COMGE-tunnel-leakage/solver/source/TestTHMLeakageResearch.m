classdef TestTHMLeakageResearch < matlab.unittest.TestCase
    % Actual public-core field/tangent/continuation and seasonal contracts.
    properties (Access=private)
        Folder
        Options
        Baseline
        FieldResult
        FieldFDResult
        ConstantResult
        DoubledKResult
    end
    methods (TestClassSetup)
        function solveSmallCoupledSystems(testCase)
            root=fileparts(mfilename('fullpath'));folder=fullfile(root,'leakage_research_20261001');
            if ~isfolder(folder),mkdir(folder);end
            testCase.Folder=tempname(folder);mkdir(testCase.Folder);
            o=joint_contact_thm_model(struct('mesh',make_thm_verification_mesh(), ...
                'durationDays',.04,'dtInitial',.02*86400,'dtMax',.02*86400, ...
                'outputEveryDays',.04,'outputPrefix',fullfile(testCase.Folder,'base'), ...
                'saveOutput',true,'verbose',false,'plot',false,'waterTableDepth',-2, ...
                'topSeepage',false,'geothermalFlux',0,'topMin',.2,'topMax',.2, ...
                'tunnelMin',-.1,'tunnelMax',-.1, ...
                'initial',struct('surfaceT',-.03,'gradient',.01), ...
                'mechanics',struct('enable',false),'solver',struct('predictor',false,'jacobianReuse',0)));
            testCase.Options=o;testCase.Baseline=TH_2D_segment_v16(o);
            o.saveOutput=false;S=o.mesh;xy=S.pAll(unique(S.tSoil(:)),:);
            o.soilKsatField=struct('coordinates_m',xy,'logMultiplier',.2*xy(:,1)-.1*xy(:,2));
            o.saveOutput=true;o.outputPrefix=fullfile(testCase.Folder,'field');
            testCase.FieldResult=TH_2D_segment_v16(o);
            o.saveOutput=false;
            o.solver.jacobian='coloredFD';testCase.FieldFDResult=TH_2D_segment_v16(o);
            o=testCase.Options;o.saveOutput=false;
            o.soilKsatField=struct('coordinates_m',[],'logMultiplier',log(2));
            testCase.ConstantResult=TH_2D_segment_v16(o);
            o=rmfield(o,'soilKsatField');o.Ksat=2*o.Ksat;
            testCase.DoubledKResult=TH_2D_segment_v16(o);
        end
    end
    methods (Test)
        function constantSpatialFieldRecoversTheUniformMaterialLaw(testCase)
            testCase.verifyEqual(testCase.ConstantResult.finalX,testCase.DoubledKResult.finalX,'AbsTol',1e-9);
            testCase.verifyEqual(testCase.ConstantResult.leakage.totalVolume_m3_per_m, ...
                testCase.DoubledKResult.leakage.totalVolume_m3_per_m,'AbsTol',1e-14);
        end
        function varyingSpatialFieldChangesTheSolvedHydraulicState(testCase)
            testCase.verifyGreaterThan(norm(testCase.FieldResult.h-testCase.Baseline.h,inf),1e-7);
            testCase.verifyGreaterThan(min(testCase.FieldResult.soilConductivityField.quadratureMultiplier),0);
            testCase.verifyEqual(testCase.FieldResult.leakage.bodyVolume_m3_per_m,0,'AbsTol',0);
        end
        function varyingFieldAnalyticTangentAgreesWithIndependentColoredFD(testCase)
            testCase.verifyEqual(testCase.FieldResult.finalX,testCase.FieldFDResult.finalX,'AbsTol',2e-6);
            testCase.verifyLessThan(max(abs(testCase.FieldResult.history(:,8))),1e-7);
            testCase.verifyLessThan(max(abs(testCase.FieldResult.history(:,9))),10);
        end
        function spatialRealizationAndInstallationSurviveARealCheckpoint(testCase)
            o=testCase.FieldResult.parameters;o.saveOutput=false;
            o.initial=struct('mode','checkpoint','checkpointFile',fullfile(testCase.Folder,'field_checkpoint.mat'));
            r=TH_2D_segment_v16(o);
            testCase.verifyEqual(r.initialX,testCase.FieldResult.finalX,'AbsTol',0);
            testCase.verifyEqual(r.soilConductivityField,testCase.FieldResult.soilConductivityField);
            testCase.verifyEqual(r.mechanics.referenceU,testCase.FieldResult.mechanics.referenceU,'AbsTol',0);
        end
        function reintegrationKeepsFullStateAndInstalledReference(testCase)
            o=testCase.continuation();o.leakage.jointFreeze.quadraturePoints=16;
            o.initial.allowNumericalReintegration=true;r=TH_2D_segment_v16(o);
            testCase.verifyEqual(r.initialX,testCase.Baseline.finalX,'AbsTol',0);
            testCase.verifyEqual(r.mechanics.referenceU,testCase.Baseline.mechanics.referenceU,'AbsTol',0);
            testCase.verifyEqual(r.mechanics.referenceTemperature,testCase.Baseline.mechanics.referenceTemperature,'AbsTol',0);
            testCase.verifyTrue(r.initialInfo.numericalReintegration.used);
        end
        function numericalFlagCannotAuthorizeAChangedSoilLaw(testCase)
            o=testCase.continuation();o.Ksat=2*o.Ksat;o.initial.allowNumericalReintegration=true;
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function numericalFlagCannotAuthorizeAChangedSpatialRealization(testCase)
            o=testCase.continuation();o.soilKsatField=struct('coordinates_m',[],'logMultiplier',.1);
            o.initial.allowNumericalReintegration=true;
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function numericalFlagCannotAuthorizeAChangedFreezingClosure(testCase)
            o=testCase.continuation();o.leakage.jointFreeze.bulkFreezingTemperature_C=-.4;
            o.initial.allowNumericalReintegration=true;
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function reintegrationCannotReuseAnOldAcceptedHistory(testCase)
            o=testCase.continuation();o.initial.allowNumericalReintegration=true;o.initial.resumeHistory=true;
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ReintegrationHistory');
        end
        function changedQuadratureStillRequiresAnExplicitNewRun(testCase)
            o=testCase.continuation();o.leakage.jointFreeze.quadraturePoints=16;
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function spatialFieldExtrapolationIsRejected(testCase)
            o=testCase.Options;o.soilKsatField=struct('coordinates_m',[0 0;1 0;0 1], ...
                'logMultiplier',[0;0;0]);
            testCase.verifyError(@()TH_2D_segment_v16(o),'SoilField:Coverage');
        end
        function seasonalBinsPreserveAnnualJointOnlyVolume(testCase)
            r=thm_seasonal_leakage_metrics(TestTHMLeakageResearch.syntheticYear());
            testCase.verifyEqual(r.annual_L_m,1e-9*365*86400*1000,'AbsTol',1e-12);
            testCase.verifyEqual(r.activeDays,365);testCase.verifyEqual(r.firstActiveDay,1);
            testCase.verifyEqual(r.lastActiveDay,365);
            testCase.verifyEqual(sum(r.monthly_L_m),r.annual_L_m,'AbsTol',1e-12);
            testCase.verifyEqual(r.signedSupplyDeficitFraction,.5,'AbsTol',1e-13);
        end
        function incompleteSeasonalHistoryIsNotAnAnnualResult(testCase)
            o=TestTHMLeakageResearch.syntheticYear();o.history(2,1)=364;o.history(2,2)=264;
            testCase.verifyError(@()thm_seasonal_leakage_metrics(o),'SeasonalLeakage:History');
        end
        function syntheticObservationsAreNeverClaimedAsExperiments(testCase)
            data=table([2;3],[98100;196200],[.1728;.3456], ...
                'VariableNames',{'temperature_C','pressureDrop_Pa','measured_L_day_m'});
            law=struct('referenceConductance_m_s',2e-10,'entryHead_m',0, ...
                'freezeSettings',joint_freeze_defaults(struct('enabled',true,'iceEnabled',false,'viscosityEnabled',false)));
            r=evaluate_joint_microchannel_observations(data,law,struct('kind','synthetic','source','unit-test'));
            testCase.verifyLessThan(r.rmse_L_day_m,1e-15);testCase.verifyFalse(r.isExperimentalComparison);
        end
        function missingObservationProvenanceIsRejected(testCase)
            data=table(2,98100,.1728,'VariableNames',{'temperature_C','pressureDrop_Pa','measured_L_day_m'});
            testCase.verifyError(@()evaluate_joint_microchannel_observations(data,struct(),struct()), ...
                'JointObservation:Provenance');
        end
        function independentCapillaryAndSlitIntegralsAgree(testCase)
            r=run_joint_leakage_reference_benchmark(fullfile(testCase.Folder,'joint_reference'));
            testCase.verifyTrue(r.passed);testCase.verifyLessThan(max(r.capillary.relativeDifference),1e-7);
            testCase.verifyLessThan(max(r.slit.relativeDifference),1e-7);
        end
    end
    methods (Access=private)
        function o=continuation(testCase)
            o=testCase.Options;o.saveOutput=false;
            o.initial=struct('mode','checkpoint','checkpointFile',[testCase.Options.outputPrefix '_checkpoint.mat']);
        end
    end
    methods (Static,Access=private)
        function o=syntheticYear()
            o=struct('history',[100 100 1e-9;365 265 1e-9], ...
                'historyColumns',{{'day','dt_day','leak_m2s'}}, ...
                'leakage',struct('jointRate_m2_per_s',[1e-9;1e-9], ...
                'bodyRate_m2_per_s',[0;0],'totalVolume_m3_per_m',1e-9*365*86400, ...
                'hydraulics',struct('unrestrictedRate_m2_per_s',[2e-9;2e-9])), ...
                'solverWork',struct('acceptedSteps',2,'rejectedSteps',0));
        end
    end
end
