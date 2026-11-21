classdef TestLocalFreezeThawDamage < matlab.unittest.TestCase
    % Literature kernels, irreversible local memory, conservative THM link.
    properties (Access=private)
        Options
        Base
        Cyclic
        NoCache
        Fine
        Folder
    end
    methods (TestClassSetup)
        function solveActualPublicCore(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.Folder=tempname(fullfile(root,'literature_damage_20261002'));
            mkdir(testCase.Folder);
            o=joint_contact_thm_model(struct('mesh',make_thm_verification_mesh(), ...
                'durationDays',30,'yearDays',10,'dtInitial',.1*86400,'dtMax',.5*86400, ...
                'outputEveryDays',5,'saveOutput',true,'plot',false,'verbose',false, ...
                'outputPrefix',fullfile(testCase.Folder,'base'),'waterTableDepth',-2, ...
                'topSeepage',false,'geothermalFlux',0,'topMin',-1,'topMax',1, ...
                'tunnelMin',-1,'tunnelMax',1,'topHTC',Inf,'tunnelHTC',Inf, ...
                'initial',struct('surfaceT',1,'gradient',0),'mechanics',struct('enable',false), ...
                'solver',struct('predictor',false,'jacobianReuse',0)));
            testCase.Base=TH_2D_segment_v16(o);
            o.leakage.jointCycleDamage=struct('enabled',true,'doseRate',.2);
            o.outputPrefix=fullfile(testCase.Folder,'cyclic');
            testCase.Options=o;testCase.Cyclic=TH_2D_segment_v16(o);
            o.saveOutput=false;o.solver.reuseConstitutiveState=false;
            testCase.NoCache=TH_2D_segment_v16(o);
            o=testCase.Options;o.saveOutput=false;o.dtMax=.25*86400;
            testCase.Fine=TH_2D_segment_v16(o);
        end
    end
    methods (Test)
        function publishedSNUsesTableCoefficientsWithoutRetuning(testCase)
            [inc,N]=luo2026_fatigue_increment([.95 .90]);
            testCase.verifyEqual(N,[6.914755724370463 65.84726351480109],'RelTol',2e-14);
            testCase.verifyEqual(inc.*N,[.05 .1],'AbsTol',1e-15);
        end
        function cyclicDamageIsNotBareMinerFailureIndex(testCase)
            [a,life]=luo2026_fatigue_increment(.95);
            testCase.verifyEqual(a,(1-.95)/life,'AbsTol',0);
            testCase.verifyNotEqual(a,1/life);
        end
        function publishedBilinearTractionAndHistory(testCase)
            m=struct('strength_Pa',10,'crackingGap_m',1,'debondingGap_m',5);
            [t,k,s]=luo2026_cohesive_mode1(3,.2,struct(),m);
            testCase.verifyEqual(t,4,'AbsTol',1e-14);testCase.verifyEqual(k,-2,'AbsTol',1e-14);
            [u,ku]=luo2026_cohesive_mode1(1,.2,s,m);
            testCase.verifyEqual(u,4/3,'AbsTol',1e-14);testCase.verifyEqual(ku,4/3,'AbsTol',1e-14);
            [v,kv]=luo2026_cohesive_mode1(-1,.2,s,m);
            testCase.verifyEqual(v,-10);testCase.verifyEqual(kv,10);
        end
        function noCycleUnderPersistentlyWarmOrColdConditions(testCase)
            c=joint_cycle_damage_defaults(struct('enabled',true));
            s=joint_freeze_thaw_memory([],0,0,true,c);
            for t=1:10,s=joint_freeze_thaw_memory(s,0,t,true,c);end
            testCase.verifyEqual(s.cycleCount,0);
            s=joint_freeze_thaw_memory([],1,0,true,c);
            s=joint_freeze_thaw_memory(s,0,1,true,c);
            testCase.verifyEqual(s.cycleCount,0); % unobserved initial freeze
        end
        function completedExcursionHasAmplitudeDependentDose(testCase)
            c=joint_cycle_damage_defaults(struct('enabled',true,'doseRate',.2));
            s=joint_freeze_thaw_memory([],[0 0],0,true,c);
            s=joint_freeze_thaw_memory(s,[.4 .8],1,true,c);
            [s,e]=joint_freeze_thaw_memory(s,[0 0],2,true,c);
            expected=([.4 .8]-.02)/.98;
            testCase.verifyEqual(e.dose,expected,'AbsTol',1e-15);
            testCase.verifyEqual(s.sampleSeverity,-expm1(-.2*expected),'AbsTol',1e-15);
            testCase.verifyEqual(s.cycleCount,[1 1]);
        end
        function chatterAndIncompleteCyclesDoNotAccumulate(testCase)
            c=joint_cycle_damage_defaults(struct('enabled',true));
            s=joint_freeze_thaw_memory([],0,0,true,c);
            for t=1:20,s=joint_freeze_thaw_memory(s,.1*(mod(t,2)),t,true,c);end
            testCase.verifyEqual(s.cumulativeDose,0);
            s=joint_freeze_thaw_memory(s,.8,21,true,c);
            s=joint_freeze_thaw_memory(s,.5,22,true,c);
            testCase.verifyTrue(s.frozen);testCase.verifyEqual(s.cycleCount,0);
        end
        function minimumDurationFiltersShortExcursions(testCase)
            c=joint_cycle_damage_defaults(struct('enabled',true,'minimumFrozenDays',2));
            s=joint_freeze_thaw_memory([],0,0,true,c);
            s=joint_freeze_thaw_memory(s,1,1,true,c);
            s=joint_freeze_thaw_memory(s,0,2,true,c);
            testCase.verifyEqual(s.cycleCount,0);testCase.verifyTrue(s.armed);
        end
        function vectorRatesAreReadyForJointVariability(testCase)
            c=joint_cycle_damage_defaults(struct('enabled',true,'doseRate',[.1;.3]));
            s=joint_freeze_thaw_memory([],zeros(2,2),0,[true;true],c);
            s=joint_freeze_thaw_memory(s,ones(2,2),1,[true;true],c);
            s=joint_freeze_thaw_memory(s,zeros(2,2),2,[true;true],c);
            testCase.verifyEqual(s.targetSeverity,-expm1(-[.1;.3]),'AbsTol',1e-15);
        end
        function coarseningConservesVolumeButChangesPhaseGeometry(testCase)
            p=testCase.Base.parameters;L=p.leakage;c=joint_cycle_damage_defaults;
            [law,info]=joint_cycle_channel_law(L,1,c);
            [~,~,i0]=joint_channel_inventory(-.1*ones(8,1),L.jointFreeze,L.jointStorage, ...
                L.jointConductance,p.liningThickness,p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
            [~,~,i1]=joint_channel_inventory(-.1*ones(8,1),law.jointFreeze,law.jointStorage, ...
                law.jointConductance,p.liningThickness,p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
            testCase.verifyEqual(i1.waterVolume_m2,i0.waterVolume_m2,'AbsTol',1e-17);
            testCase.verifyLessThan(i1.liquidVolumeFraction,i0.liquidVolumeFraction);
            testCase.verifyEqual(info.channelDensityMultiplier,ones(8,1)/2.25,'AbsTol',1e-15);
        end
        function actualTHMHasDamageAndMoreThawedLeakage(testCase)
            r=testCase.Cyclic;
            testCase.verifyGreaterThan(r.jointCycleDamage.appliedSeverity(7),.1);
            testCase.verifyGreaterThan(r.leakage.totalVolume_m3_per_m, ...
                testCase.Base.leakage.totalVolume_m3_per_m);
            testCase.verifyGreaterThan(norm(r.finalX-testCase.Base.finalX,inf),1e-8);
            testCase.verifyEqual(r.leakage.bodyVolume_m3_per_m,0);
            testCase.verifyEqual(r.leakage.volumeByJoint_m3_per_m([1:6 8]),zeros(1,7));
        end
        function acceptedWaterAndHeatBudgetsRemainStrict(testCase)
            r=testCase.Cyclic;
            testCase.verifyLessThan(max(abs(r.history(:,8))),1e-7);
            testCase.verifyLessThan(max(abs(r.history(:,9))),10);
            testCase.verifyTrue(audit_joint_cycle_damage(r).passed);
        end
        function disablingStateCacheDoesNotEraseOldGeometry(testCase)
            testCase.verifyEqual(testCase.NoCache.finalX,testCase.Cyclic.finalX,'AbsTol',1e-9);
            testCase.verifyEqual(testCase.NoCache.jointCycleDamage.history, ...
                testCase.Cyclic.jointCycleDamage.history,'AbsTol',1e-10);
        end
        function fineStepDoesNotLoseCycles(testCase)
            c=testCase.Cyclic;f=testCase.Fine;
            testCase.verifyEqual(f.jointCycleDamage.finalMemory.cycleCount, ...
                c.jointCycleDamage.finalMemory.cycleCount);
            testCase.verifyLessThan(abs(f.leakage.totalVolume_m3_per_m/c.leakage.totalVolume_m3_per_m-1),.05);
        end
        function disabledEvolutionPreservesThePriorPhysicalTrajectory(testCase)
            o=testCase.Options;o.saveOutput=false;o.leakage.jointCycleDamage.enabled=false;
            r=TH_2D_segment_v16(o);
            testCase.verifyEqual(r.finalX,testCase.Base.finalX,'AbsTol',0);
            testCase.verifyEqual(r.history,testCase.Base.history,'AbsTol',0);
            old=testCase.Base.parameters;old.leakage=rmfield(old.leakage,'jointCycleDamage');
            testCase.verifyEqual(thm_constitutive_signature(old), ...
                thm_constitutive_signature(testCase.Base.parameters));
        end
        function completeMemorySurvivesConditionalCheckpoint(testCase)
            o=testCase.Options;o.saveOutput=false;o.durationDays=1;
            o.initial=struct('mode','checkpoint','checkpointFile', ...
                fullfile(testCase.Folder,'cyclic_checkpoint.mat'));
            r=TH_2D_segment_v16(o);
            testCase.verifyEqual(r.initialX,testCase.Cyclic.finalX,'AbsTol',0);
            testCase.verifyEqual(r.jointCycleDamage.initialMemory, ...
                testCase.Cyclic.jointCycleDamage.finalMemory);
            testCase.verifyEqual(r.mechanics.referenceU,testCase.Cyclic.mechanics.referenceU);
        end
        function alteredDamageLawCannotReuseOldCheckpoint(testCase)
            o=testCase.Options;o.saveOutput=false;o.leakage.jointCycleDamage.doseRate=.3;
            o.initial=struct('mode','checkpoint','checkpointFile', ...
                fullfile(testCase.Folder,'cyclic_checkpoint.mat'),'allowNumericalReintegration',true);
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function corruptedInternalDamageIsRejected(testCase)
            s=testCase.Cyclic.jointCycleDamage.finalMemory;s.cumulativeDose(7,1)=s.cumulativeDose(7,1)+1;
            testCase.verifyError(@()validate_joint_cycle_memory(s,8,size(s.lastIce,2), ...
                testCase.Cyclic.jointCycleDamage.settings),'JointCycle:Memory');
        end
        function partialFrozenExcursionSurvivesCheckpoint(testCase)
            o=testCase.Options;o.dtInitial=.5*86400;o.outputPrefix=fullfile(testCase.Folder,'split');
            o.durationDays=1;first=TH_2D_segment_v16(o);
            testCase.verifyTrue(any(first.jointCycleDamage.finalMemory.frozen,'all'));
            o.initial=struct('mode','checkpoint','checkpointFile',[o.outputPrefix '_checkpoint.mat']);
            o.durationDays=29;o.saveOutput=false;last=TH_2D_segment_v16(o);
            o=testCase.Options;o.dtInitial=.5*86400;o.saveOutput=false;full=TH_2D_segment_v16(o);
            testCase.verifyEqual(last.finalX,full.finalX,'AbsTol',1e-10);
            testCase.verifyEqual(last.jointCycleDamage.finalMemory,full.jointCycleDamage.finalMemory);
            testCase.verifyEqual(last.mechanics.referenceU,full.mechanics.referenceU);
        end
        function analyticTangentIncludesChangingPhaseDistribution(testCase)
            o=testCase.Options;o.durationDays=8;o.saveOutput=false;
            a=TH_2D_segment_v16(o);o.solver.jacobian='coloredFD';f=TH_2D_segment_v16(o);
            testCase.verifyGreaterThan(a.jointCycleDamage.appliedSeverity(7),0);
            testCase.verifyEqual(a.finalX,f.finalX,'AbsTol',2e-6);
            testCase.verifyEqual(a.jointCycleDamage.history(:,2:17), ...
                f.jointCycleDamage.history(:,2:17),'AbsTol',1e-6);
        end
        function scalarConductanceHasTheSameEightJointTrajectory(testCase)
            o=testCase.Options;o.durationDays=.1;o.saveOutput=false;
            % Isolate existing-channel input shape; no additional open slit.
            o.thm.openingTransmissionFraction=0;
            o.leakage.jointConductance=1e-12;
            a=TH_2D_segment_v16(o);o.leakage.jointConductance=ones(8,1)*1e-12;
            b=TH_2D_segment_v16(o);
            testCase.verifyEqual(a.finalX,b.finalX,'AbsTol',0);
            testCase.verifyEqual(a.jointCycleDamage.history,b.jointCycleDamage.history,'AbsTol',0);
        end
    end
end
