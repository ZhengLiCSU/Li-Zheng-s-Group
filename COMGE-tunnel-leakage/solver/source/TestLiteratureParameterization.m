classdef TestLiteratureParameterization < matlab.unittest.TestCase
    % Public unit mapping, fitting identifiability and internal-state guards.
    methods (Test)
        function kanagawaParametersUsePublishedUnits(testCase)
            [p,e]=thm_literature_parameter_pack('kanagawa_th');
            testCase.verifyEqual([p.porosity,p.thetaR,p.vgAlpha,p.vgN,p.Ksat,p.ks], ...
                [.535,.05,1.11,1.48,3.2e-6,2.5],'AbsTol',1e-14);
            testCase.verifyEqual(p.rhoS*p.cs,2.3e6,'AbsTol',1e-9);
            testCase.verifyEqual(height(e.ledger),7);
            testCase.verifyFalse(e.isSiteCalibration);
            testCase.verifyFalse(e.isJointDamageCalibration);
        end
        function siltTupleDoesNotMixPorosityAndRetention(testCase)
            [p,e]=thm_literature_parameter_pack('silt_retention');
            testCase.verifyEqual([p.porosity,p.thetaR,p.vgAlpha,p.vgN], ...
                [.49,.05,.65,1.67],'AbsTol',1e-15);
            testCase.verifyFalse(isfield(p,'Ksat'));
            testCase.verifyTrue(any(contains(string(e.missingCalibration),'Ksat')));
        end
        function kanagawaWarmInitialRetentionMatchesPublishedValue(testCase)
            p=thm_literature_parameter_pack();m=1-1/p.vgN;
            theta=p.thetaR+(p.porosity-p.thetaR)*(1+(p.vgAlpha*2.46675)^p.vgN)^(-m);
            testCase.verifyEqual(theta,.33,'AbsTol',1e-6);
        end
        function unsupportedSoilPackIsRejected(testCase)
            testCase.verifyError(@()thm_literature_parameter_pack('unmeasured_joint'), ...
                'THMLiterature:Name');
        end
        function syntheticParametersRecoveredWithBroadDoseSpan(testCase)
            [d,p]=TestLiteratureParameterization.syntheticData([0;1;2;5;10;20;40],.05,.5);
            r=fit_joint_cycle_coarsening(d,p);
            testCase.verifyTrue(r.identifiable);
            testCase.verifyEqual(r.doseRate,.05,'RelTol',1e-6);
            testCase.verifyEqual(r.radiusGrowthAtUnitDamage,.5,'RelTol',1e-6);
            testCase.verifyTrue(r.isSynthetic);
            testCase.verifyFalse(r.isIndependentValidation);
        end
        function narrowEarlyDoseSpanDoesNotIdentifyTwoParameters(testCase)
            [d,p]=TestLiteratureParameterization.syntheticData([0;1e-4;2e-4;5e-4;1e-3],.05,.5);
            r=fit_joint_cycle_coarsening(d,p);
            testCase.verifyFalse(r.identifiable);
            testCase.verifyNotEmpty(r.warnings);
        end
        function concreteReferenceCannotBecomeGasketParameterEvidence(testCase)
            [d,p]=TestLiteratureParameterization.syntheticData([0;1;2;5;10;20;40],.05,.5);
            p.material='concrete_joint';p.kind='literature_table';
            r=fit_joint_cycle_coarsening(d,p);
            testCase.verifyFalse(r.materialMatchesSegmentGasket);
            testCase.verifyTrue(any(contains(string(r.warnings),'must not be imported')));
        end
        function changingCompressionAndUnthawedObservationsRejected(testCase)
            [d,p]=TestLiteratureParameterization.syntheticData([0;1;2;5;10;20;40],.05,.5);
            p.constantCompression=false;
            testCase.verifyError(@()fit_joint_cycle_coarsening(d,p),'JointCycleFit:Conditions');
            p.constantCompression=true;p.thawedReferenceCorrected=false;
            testCase.verifyError(@()fit_joint_cycle_coarsening(d,p),'JointCycleFit:Conditions');
        end
        function zeroDoseReferenceIsRequired(testCase)
            [d,p]=TestLiteratureParameterization.syntheticData([1;2;5;10;20;40],.05,.5);
            testCase.verifyError(@()fit_joint_cycle_coarsening(d,p),'JointCycleFit:Coverage');
        end
        function heterogeneousMeanDoseCannotReplaceLocalPathExposure(testCase)
            [d,p]=TestLiteratureParameterization.syntheticData([0;1;2;5;10;20;40],.05,.5);
            p.uniformPathDose=false;
            testCase.verifyError(@()fit_joint_cycle_coarsening(d,p),'JointCycleFit:Conditions');
        end
        function damageLawHasStableSmallDoseAndCorrectDerivatives(testCase)
            [s,k,x]=joint_cycle_severity([1e-12 2],[.2 .2],.1);
            testCase.verifyEqual(s(1)-.1,.9*.2e-12,'AbsTol',1e-17);
            testCase.verifyEqual(k,.9*[1e-12 2].*exp(-.2*[1e-12 2]),'AbsTol',1e-15);
            testCase.verifyEqual(x,.9*.2*exp(-.2*[1e-12 2]),'AbsTol',1e-15);
        end
        function invalidActiveMaskCannotSilentlyActivateASealedJoint(testCase)
            c=joint_cycle_damage_defaults();
            testCase.verifyError(@()joint_freeze_thaw_memory([],0,0,NaN,c),'JointCycle:Mask');
            testCase.verifyError(@()joint_freeze_thaw_memory([],0,0,2,c),'JointCycle:Mask');
        end
        function impossibleDoseCannotBeHiddenByConsistentDamage(testCase)
            c=joint_cycle_damage_defaults();s=joint_freeze_thaw_memory([],0,0,true,c);
            s.cumulativeDose=2;s.sampleSeverity=joint_cycle_severity(2,c.doseRate,0);
            s.targetSeverity=s.sampleSeverity;
            testCase.verifyError(@()validate_joint_cycle_memory(s,1,1,c),'JointCycle:Memory');
        end
        function corruptedFrozenPeakIsRejected(testCase)
            c=joint_cycle_damage_defaults();s=joint_freeze_thaw_memory([],0,0,true,c);
            s=joint_freeze_thaw_memory(s,.8,1,true,c);s.peakIce=.5;
            testCase.verifyError(@()joint_freeze_thaw_memory(s,.9,2,true,c),'JointCycle:Memory');
        end
    end
    methods (Static,Access=private)
        function [d,p]=syntheticData(dose,rate,alpha)
            ratio=(1+alpha*(-expm1(-rate*dose))).^2;
            d=table(dose,ratio,'VariableNames',{'cumulativeDose','conductanceRatio'});
            p=struct('kind','synthetic','source','analytic recovery test only', ...
                'material','gasketed_segment_joint','thawedReferenceCorrected',true, ...
                'constantCompression',true,'uniformPathDose',true);
        end
    end
end
