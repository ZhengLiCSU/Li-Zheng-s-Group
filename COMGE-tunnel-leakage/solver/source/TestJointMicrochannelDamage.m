classdef TestJointMicrochannelDamage < matlab.unittest.TestCase
    % Physical moment, phase, coupled tangent and continuation contracts.
    properties (Access=private)
        Folder
        Options
        Base
        Zero
        Damaged
        Explicit
        FiniteDifference
    end
    methods (TestClassSetup)
        function solvePublicCoreCases(testCase)
            root=fileparts(mfilename('fullpath'));
            folder=fullfile(root,'joint_damage_20261002');
            if ~isfolder(folder),mkdir(folder);end
            testCase.Folder=tempname(folder);mkdir(testCase.Folder);
            o=joint_contact_thm_model(struct('mesh',make_thm_verification_mesh(), ...
                'durationDays',.04,'dtInitial',.02*86400,'dtMax',.02*86400, ...
                'outputEveryDays',.04,'saveOutput',true,'plot',false,'verbose',false, ...
                'outputPrefix',fullfile(testCase.Folder,'base'), ...
                'waterTableDepth',-2,'topSeepage',false,'geothermalFlux',0, ...
                'topMin',.2,'topMax',.2,'tunnelMin',-.1,'tunnelMax',-.1, ...
                'initial',struct('surfaceT',-.03,'gradient',.01), ...
                'mechanics',struct('enable',false), ...
                'solver',struct('predictor',false,'jacobianReuse',0)));
            testCase.Options=o;testCase.Base=TH_2D_segment_v16(o);
            o.leakage.jointDamage=struct('enabled',true,'severity',0);
            o.saveOutput=false;testCase.Zero=TH_2D_segment_v16(o);
            o.leakage.jointDamage.severity=.5;
            o.saveOutput=true;o.outputPrefix=fullfile(testCase.Folder,'damaged');
            testCase.Damaged=TH_2D_segment_v16(o);
            [law,~]=joint_microchannel_damage(testCase.Damaged.parameters.leakage);
            o=testCase.Options;o.saveOutput=false;o.leakage=law;
            o.leakage.jointDamage.enabled=false;
            testCase.Explicit=TH_2D_segment_v16(o);
            o=testCase.Damaged.parameters;o.saveOutput=false;o.solver.jacobian='coloredFD';
            testCase.FiniteDifference=TH_2D_segment_v16(o);
        end
    end
    methods (Test)
        function zeroDamagePreservesTheActualCoupledTrajectory(testCase)
            testCase.verifyEqual(testCase.Zero.finalX,testCase.Base.finalX,'AbsTol',0);
            testCase.verifyEqual(testCase.Zero.history,testCase.Base.history,'AbsTol',0);
        end
        function resolvedDamageMatchesAnExplicitPhysicalGeometry(testCase)
            testCase.verifyEqual(testCase.Damaged.finalX,testCase.Explicit.finalX,'AbsTol',0);
            testCase.verifyEqual(testCase.Damaged.history,testCase.Explicit.history,'AbsTol',0);
        end
        function analyticTangentMatchesColoredFiniteDifferences(testCase)
            testCase.verifyEqual(testCase.Damaged.finalX,testCase.FiniteDifference.finalX,'AbsTol',2e-6);
        end
        function thermalInventoryChangesWithRadiusSquared(testCase)
            lambda=1.25;
            testCase.verifyEqual(testCase.Damaged.state.jointInventory.waterVolume_m2, ...
                testCase.Base.state.jointInventory.waterVolume_m2*lambda^2,'RelTol',1e-13);
            testCase.verifyEqual(testCase.Damaged.jointDamage.referenceConductanceMultiplier, ...
                lambda^4*ones(8,1),'AbsTol',1e-14);
        end
        function damageChangesTheSolvedStateRatherThanScalingOutput(testCase)
            testCase.verifyGreaterThan(norm(testCase.Damaged.finalX-testCase.Base.finalX,inf),1e-7);
            testCase.verifyLessThan(max(abs(testCase.Damaged.history(:,8))),1e-7);
            testCase.verifyLessThan(max(abs(testCase.Damaged.history(:,9))),10);
            testCase.verifyEqual(testCase.Damaged.leakage.bodyVolume_m3_per_m,0,'AbsTol',0);
        end
        function parameterReplayDoesNotApplyGeometryTwice(testCase)
            o=testCase.Damaged.parameters;o.saveOutput=false;
            replay=TH_2D_segment_v16(o);
            testCase.verifyEqual(replay.finalX,testCase.Damaged.finalX,'AbsTol',0);
        end
        function physicalDamageAndInstallationSurviveCheckpoint(testCase)
            o=testCase.Damaged.parameters;o.saveOutput=false;
            o.initial=struct('mode','checkpoint','checkpointFile',fullfile(testCase.Folder,'damaged_checkpoint.mat'));
            r=TH_2D_segment_v16(o);
            testCase.verifyEqual(r.initialX,testCase.Damaged.finalX,'AbsTol',0);
            testCase.verifyEqual(r.jointDamage,testCase.Damaged.jointDamage);
            testCase.verifyEqual(r.mechanics.referenceU,testCase.Damaged.mechanics.referenceU,'AbsTol',0);
        end
        function changedDamageCannotReusePhysicalHistory(testCase)
            o=testCase.Damaged.parameters;o.saveOutput=false;
            o.leakage.jointDamage.severity=.7;
            o.initial=struct('mode','checkpoint','checkpointFile',fullfile(testCase.Folder,'damaged_checkpoint.mat'), ...
                'allowNumericalReintegration',true);
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function disabledDamageRetainsLegacySignature(testCase)
            p=testCase.Base.parameters;
            old=p;old.leakage=rmfield(old.leakage,'jointDamage');
            testCase.verifyEqual(thm_constitutive_signature(p),thm_constitutive_signature(old));
        end
        function aSealedJointDoesNotAcquireInventedConnectivity(testCase)
            L=testCase.Base.parameters.leakage;
            L.jointDamage=struct('enabled',true,'severity',1);
            law=joint_microchannel_damage(L);
            testCase.verifyEqual(law.jointConductance(L.jointConductance==0), ...
                zeros(nnz(L.jointConductance==0),1),'AbsTol',0);
        end
        function largerPoresFreezeMoreReadily(testCase)
            L=testCase.Base.parameters.leakage;
            L.jointDamage=struct('enabled',true,'severity',1);
            law=joint_microchannel_damage(L);
            f0=joint_freeze_conductance(-.1*ones(8,1),L.jointFreeze);
            f1=joint_freeze_conductance(-.1*ones(8,1),law.jointFreeze);
            testCase.verifyLessThan(f1,f0);
        end
        function deeplyFrozenLeakageIsNotForcedToIncrease(testCase)
            L=testCase.Base.parameters.leakage;
            L.jointDamage=struct('enabled',true,'severity',1);
            law=joint_microchannel_damage(L);
            c0=L.jointConductance(:).*joint_freeze_conductance(-2*ones(8,1),L.jointFreeze);
            c1=law.jointConductance(:).*joint_freeze_conductance(-2*ones(8,1),law.jointFreeze);
            testCase.verifyLessThan(c1(7),c0(7));
        end
        function channelDensityGrowthChangesBothFlowAndStorage(testCase)
            L=testCase.Base.parameters.leakage;
            L.jointDamage=struct('enabled',true,'severity',.5, ...
                'radiusGrowthAtUnitDamage',0,'channelDensityGrowthAtUnitDamage',2);
            [law,info]=joint_microchannel_damage(L);
            testCase.verifyEqual(law.jointConductance,2*L.jointConductance(:),'AbsTol',1e-20);
            testCase.verifyEqual(info.connectedVolumeMultiplier,2*ones(8,1),'AbsTol',0);
            testCase.verifyEqual(law.jointFreeze.flowWeightedRadius_m,L.jointFreeze.flowWeightedRadius_m,'AbsTol',0);
        end
        function invalidDamageCannotBeClipped(testCase)
            L=testCase.Base.parameters.leakage;L.jointDamage=struct('enabled',true,'severity',1.1);
            testCase.verifyError(@()joint_microchannel_damage(L),'JointDamage:Range');
        end
        function missingInventoryIsRejected(testCase)
            L=testCase.Base.parameters.leakage;L.jointStorage.enabled=false;
            L.jointDamage=struct('enabled',true,'severity',.5);
            testCase.verifyError(@()joint_microchannel_damage(L),'JointDamage:Inventory');
        end
        function enlargedChannelsStillObeyThePhysicalVolumeCeiling(testCase)
            o=testCase.Options;o.saveOutput=false;
            o.leakage.jointDamage=struct('enabled',true,'severity',1,'radiusGrowthAtUnitDamage',1e6);
            testCase.verifyError(@()TH_2D_segment_v16(o),'TH:channelGeometry');
        end
    end
end
