classdef TestTHMJointFreezingContinuation < matlab.unittest.TestCase
    % Actual checkpoint compatibility for the corrected nonzero joint point.
    properties (Access=private)
        Folder
        Options
        First
        Checkpoint
    end
    methods (TestClassSetup)
        function createOneCorrectedCheckpoint(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.Folder=tempname(fullfile(root,'physics_validation_20261001'));
            mkdir(testCase.Folder);
            testCase.Options=joint_supply_th_model(struct( ...
                'mesh',make_thm_verification_mesh(),'durationDays',.02, ...
                'dtInitial',.01*86400,'dtMax',.01*86400,'outputEveryDays',.02, ...
                'outputPrefix',fullfile(testCase.Folder,'corrected'), ...
                'plot',false,'verbose',false,'saveOutput',true,'waterTableDepth',-2, ...
                'topSeepage',false,'geothermalFlux',0,'topMin',1,'topMax',1, ...
                'tunnelMin',1,'tunnelMax',1,'initial',struct('surfaceT',1,'gradient',0), ...
                'mechanics',struct('enable',false), ...
                'leakage',struct('jointFreeze',struct('bulkFreezingTemperature_C',-.4))));
            testCase.First=TH_2D_segment_v16(testCase.Options);
            testCase.Checkpoint=load([testCase.Options.outputPrefix '_checkpoint.mat']);
        end
    end
    methods (Test)
        function newNonzeroPointCheckpointHasCorrectClosureGuard(testCase)
            cp=testCase.Checkpoint.checkpoint;
            testCase.verifyEqual(cp.signature.jointSlitFreezingClosureRevision, ...
                'common_bulk_freezing_point_20261001');
            testCase.verifyEqual(cp.signature.leakage.jointFreeze.bulkFreezingTemperature_C, ...
                -.4,'AbsTol',0);
        end
        function legacyNonzeroPointCheckpointCannotBeSilentlyContinued(testCase)
            data=testCase.Checkpoint;
            data.checkpoint.signature=rmfield(data.checkpoint.signature, ...
                'jointSlitFreezingClosureRevision');
            file=fullfile(testCase.Folder,'legacy_nonzero_point_checkpoint.mat');
            save(file,'-struct','data','-v7');
            o=testCase.Options;o.saveOutput=false;
            o.initial=struct('mode','checkpoint','checkpointFile',file);
            testCase.verifyError(@()TH_2D_segment_v16(o),'THM:ConstitutiveMismatch');
        end
        function correctedCheckpointContinuesWithTheOriginalInstallation(testCase)
            o=testCase.Options;o.saveOutput=false;
            o.initial=struct('mode','checkpoint', ...
                'checkpointFile',[testCase.Options.outputPrefix '_checkpoint.mat']);
            continued=TH_2D_segment_v16(o);
            testCase.verifyEqual(continued.initialX,testCase.First.finalX,'AbsTol',0);
            testCase.verifyEqual(continued.mechanics.referenceU, ...
                testCase.First.mechanics.referenceU,'AbsTol',0);
            testCase.verifyEqual(continued.mechanics.referenceTemperature, ...
                testCase.First.mechanics.referenceTemperature,'AbsTol',0);
        end
        function zeroPointKeepsTheExistingDefaultSignature(testCase)
            p=testCase.First.parameters;p.leakage.jointFreeze.bulkFreezingTemperature_C=0;
            old=rmfield(p,{'mesh','meshFile','durationDays','dtInitial','dtMax','dtMin', ...
                'outputEveryDays','outputPrefix','saveOutput','plot','verbose','initial','solver','spinup','soilKsatField'});
            old.leakage=rmfield(old.leakage,{'jointDamage','jointCycleDamage'});
            testCase.verifyEqual(thm_constitutive_signature(p),old);
        end
    end
end
