classdef TestV16THM < matlab.unittest.TestCase
    properties
        coupledReport
    end
    methods (TestClassSetup)
        function exerciseMonolithicSystem(testCase)
            testCase.coupledReport=TH_2D_segment_v16('thmselftest');
        end
    end
    methods (Test)
        function completeTangentHasAllCoupledDependencies(testCase)
            testCase.verifyLessThan(testCase.coupledReport.jacobianRelativeError,1e-5);
        end
        function deformationChangesWaterSolution(testCase)
            testCase.verifyGreaterThan(testCase.coupledReport.maximumHeadFeedback_m,1e-7);
            testCase.verifyGreaterThan(testCase.coupledReport.maxPorosityChange,1e-9);
        end
        function realOpenBranchAndTrialDomainAreExercised(testCase)
            testCase.verifyGreaterThan(testCase.coupledReport.trialOpenConductance_m_s,0);
            testCase.verifyTrue(testCase.coupledReport.oversizedTrialRejected);
            testCase.verifyLessThan(testCase.coupledReport.disabledIceViscositySwitchError,1e-13);
        end
        function movingPoresPreserveEnergyDatumInvariance(testCase)
            testCase.verifyLessThan(testCase.coupledReport.enthalpyDatumRelativeError,1e-9);
        end
        function bothJointBranchesRespectThePrescribedFreezingPoint(testCase)
            testCase.verifyEqual(testCase.coupledReport.depressedFreezingSubzeroSlitMobility, ...
                1,'AbsTol',1e-13);
            testCase.verifyEqual(testCase.coupledReport.depressedFreezingColdSlitConductance_m_s, ...
                0,'AbsTol',0);
        end
        function acceptedStepsConserveWaterHeatAndMechanicalBalance(testCase)
            testCase.verifyLessThan(testCase.coupledReport.maximumWaterDefect_m2,1e-7);
            testCase.verifyLessThan(testCase.coupledReport.maximumEnergyDefect_J_per_m,10);
            testCase.verifyLessThan(testCase.coupledReport.maximumFreeMechanicalResidual_Nm,.01);
            testCase.verifyLessThan(testCase.coupledReport.maxVolumetricStrain,.05);
            out=testCase.coupledReport.coupledOutput;
            testCase.verifyEqual(size(out.mechanics.jointHistory,1),out.solverWork.acceptedSteps);
            testCase.verifyEqual(size(out.mechanics.jointHistory,2),73);
            testCase.verifyEqual(size(out.mechanics.strainHistory),[out.solverWork.acceptedSteps,7]);
            s=out.snapshots(end);
            testCase.verifyEqual([s.T;s.h;s.U],out.finalX,'AbsTol',0);
            testCase.verifyEqual(out.leakage.bodyVolume_m3_per_m,0,'AbsTol',0);
            testCase.verifyEqual(out.leakage.volumeByJoint_m3_per_m([1:6,8]),zeros(1,7),'AbsTol',0);
        end
        function checkpointRetainsInstalledMechanicalReference(testCase)
            r=TestV16THM.restartExercise();
            testCase.verifyEqual(r.referenceU,r.resumedReferenceU,'AbsTol',0);
            testCase.verifyEqual(r.referenceT,r.resumedReferenceT,'AbsTol',0);
            testCase.verifyEqual(r.initialX,r.checkpointX,'AbsTol',0);
            testCase.verifyEqual(r.finalX,r.uninterruptedX,'AbsTol',2e-7);
            testCase.verifyEqual(r.referenceHistoryRows,r.referenceAcceptedSteps);
        end
    end
    methods (Static,Access=private)
        function r=restartExercise()
            folder=tempname(fileparts(mfilename('fullpath'))); mkdir(folder);
            o=joint_supply_th_model(struct('mesh',make_thm_verification_mesh(), ...
                'durationDays',.1,'dtInitial',.025*86400,'dtMax',.025*86400, ...
                'outputEveryDays',.05,'plot',false,'saveOutput',false,'verbose',false, ...
                'waterTableDepth',-2,'topSeepage',false,'geothermalFlux',0, ...
                'topMin',.2,'topMax',.2,'tunnelMin',-.1,'tunnelMax',-.1, ...
                'initial',struct('surfaceT',-.03,'gradient',.01), ...
                'mechanics',struct('enable',false), ...
                'solver',struct('predictor',false,'jacobianReuse',0)));
            reference=TH_2D_segment_v16(o);
            o.durationDays=.05; o.saveOutput=true; o.outputPrefix=fullfile(folder,'partial');
            first=TH_2D_segment_v16(o); file=[o.outputPrefix '_checkpoint.mat'];
            d=load(file,'checkpoint','mechanicalHistory','jointMechanicalHistory');
            o.initial=struct('mode','checkpoint','checkpointFile',file);
            o.outputPrefix=fullfile(folder,'resumed');
            second=TH_2D_segment_v16(o);
            r=struct('referenceU',first.mechanics.referenceU, ...
                'resumedReferenceU',second.mechanics.referenceU, ...
                'referenceT',first.mechanics.referenceTemperature, ...
                'resumedReferenceT',second.mechanics.referenceTemperature, ...
                'initialX',second.initialX,'checkpointX',d.checkpoint.X, ...
                'finalX',second.finalX,'uninterruptedX',reference.finalX, ...
                'referenceHistoryRows',size(d.jointMechanicalHistory,1), ...
                'referenceAcceptedSteps',first.solverWork.acceptedSteps);
        end
    end
end
