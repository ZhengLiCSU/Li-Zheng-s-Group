classdef TestTHMCompletedRecovery < matlab.unittest.TestCase
    % Complete-year recovery performs zero new steps and validates archives.
    % One small warm/sealed year supplies the ORIGINAL accepted full history.
    % Separate copies exercise the final-CP branch and the month-before-cycle
    % crash window. No real campaign checkpoint is edited or used as a fixture.
    properties (Access=private)
        Evidence
    end

    methods (TestClassSetup)
        function prepareCompletedYearEvidence(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.Evidence=TestTHMCompletedRecovery.exercise(root);
        end
    end

    methods (Test)
        function completeCoreHistoryDoesNotAddStepsOrDuplicateCycles(testCase)
            r=testCase.Evidence; actual=r.completeCore; expected=r.sourceOut;
            testCase.verifyEqual(actual.initialX,expected.initialX,'AbsTol',0);
            testCase.verifyEqual(actual.finalX,expected.finalX,'AbsTol',0);
            testCase.verifyEqual(actual.history,expected.history,'AbsTol',0);
            testCase.verifyEqual(actual.cycles,expected.cycles,'AbsTol',0);
            testCase.verifySize(actual.cycles,[1,5]);
            testCase.verifyEqual(actual.solverWork.acceptedSteps,expected.solverWork.acceptedSteps);
            testCase.verifyEqual(actual.mechanics.referenceU,expected.mechanics.referenceU,'AbsTol',0);
            testCase.verifyEqual(actual.mechanics.referenceTemperature,expected.mechanics.referenceTemperature,'AbsTol',0);
            testCase.verifyEqual(actual.timeDay,365,'AbsTol',0);
        end

        function missingYearEndResultsAndCycleAreRebuiltWithoutASecondYear(testCase)
            r=testCase.Evidence; actual=r.missingResults.lastYear; expected=r.sourceOut;
            testCase.verifyEqual(r.missingResults.completedYears,1);
            testCase.verifyEqual(actual.finalX,expected.finalX,'AbsTol',0);
            testCase.verifyEqual(actual.initialX,expected.initialX,'AbsTol',0);
            testCase.verifyEqual(actual.history,expected.history,'AbsTol',0);
            testCase.verifyEqual(actual.cycles,expected.cycles,'AbsTol',0);
            testCase.verifySize(actual.cycles,[1,5]);
            testCase.verifyEqual(actual.solverWork.acceptedSteps,expected.solverWork.acceptedSteps);
            testCase.verifyEqual(actual.mechanics.referenceU,expected.mechanics.referenceU,'AbsTol',0);
            testCase.verifyEqual(actual.mechanics.referenceTemperature,expected.mechanics.referenceTemperature,'AbsTol',0);
            testCase.verifyEqual(actual.timeDay,365,'AbsTol',0);
            testCase.verifyEqual(sum(actual.history(:,2)),365,'AbsTol',1e-10);
            testCase.verifySize(r.missingResults.annual,[1,14]);
            testCase.verifySize(r.missingResults.monthly,[12,11]);
            testCase.verifySize(r.missingResults.daily,[365,12]);
            testCase.verifySize(r.missingResults.byJoint,[8,4]);
            testCase.verifyEqual(actual.leakage.volumeByJoint_m3_per_m,zeros(1,8),'AbsTol',0);
            testCase.verifyEqual(actual.leakage.totalVolume_m3_per_m,0,'AbsTol',0);
            testCase.verifyTrue(actual.parameters.initial.resumeHistory);
            testCase.verifyTrue(isfile(fullfile(r.missingResults.config.outputFolder,'year_01_results.mat')));
        end

        function existingResultsAreVerifiedByTheCoreAndNotTrustOnlyRecovered(testCase)
            r=testCase.Evidence; actual=r.existingResults.lastYear; expected=r.sourceOut;
            testCase.verifyEqual(r.existingResults.completedYears,1);
            testCase.verifyEqual(actual.finalX,expected.finalX,'AbsTol',0);
            testCase.verifyEqual(actual.initialX,expected.initialX,'AbsTol',0);
            testCase.verifyEqual(actual.history,expected.history,'AbsTol',0);
            testCase.verifyEqual(actual.cycles,expected.cycles,'AbsTol',0);
            testCase.verifyEqual(actual.solverWork.acceptedSteps,expected.solverWork.acceptedSteps);
            testCase.verifyTrue(actual.parameters.initial.resumeHistory);
            testCase.verifyEqual(actual.mesh,expected.mesh);
            testCase.verifyEqual(actual.pAll,expected.pAll,'AbsTol',0);
            testCase.verifyEqual(actual.soilNodes,expected.soilNodes,'AbsTol',0);
            testCase.verifyEqual(actual.mechanics.history,expected.mechanics.history,'AbsTol',0);
            testCase.verifyEqual(actual.leakage.jointRate_m2_per_s,expected.leakage.jointRate_m2_per_s,'AbsTol',0);
        end

        function recoveredArchiveWithChangedPhysicsIsRejected(testCase)
            r=testCase.Evidence;
            testCase.verifyError(@() TestTHMCompletedRecovery.corruptedArchive(r,'physical'), ...
                'THMAnnual:RecoveredPhysicalConfig');
        end

        function recoveredArchiveWithChangedMeshIsRejected(testCase)
            r=testCase.Evidence;
            testCase.verifyError(@() TestTHMCompletedRecovery.corruptedArchive(r,'mesh'), ...
                'THMAnnual:RecoveredMesh');
        end

        function unreadableCurrentCheckpointUsesCompatibleAlternate(testCase)
            r=testCase.Evidence; actual=r.alternateResults.lastYear;
            testCase.verifyEqual(actual.history,r.sourceOut.history,'AbsTol',0);
            testCase.verifyEqual(actual.finalX,r.sourceOut.finalX,'AbsTol',0);
            testCase.verifyEqual(actual.initialX,r.sourceOut.initialX,'AbsTol',0);
            testCase.verifyEqual(actual.cycles,r.sourceOut.cycles,'AbsTol',0);
            testCase.verifyEqual(actual.parameters.initial.checkpointFile,r.alternateCheckpointFile);
            testCase.verifyEqual(r.alternateResults.completedYears,1);
            testCase.verifyTrue(isfile(r.alternateCheckpointFile));
        end

        function readableIncompatibleCheckpointIsStillRejected(testCase)
            r=testCase.Evidence;
            testCase.verifyError(@() TestTHMCompletedRecovery.corruptedArchive(r,'version'), ...
                'THMAnnual:LegacyFolder');
        end
    end

    methods (Static,Access=private)
        function r=exercise(root)
            parent=fullfile(root,'thm_completed_recovery_checks');
            if ~isfolder(parent), mkdir(parent); end
            folder=tempname(parent); mkdir(folder);
            sourceFolder=fullfile(folder,'source_year');
            model=struct('mesh',make_thm_verification_mesh(16), ...
                'dtInitial',4*86400,'dtMax',4*86400,'dtMin',86400/16384, ...
                'waterTableDepth',-2,'topSeepage',false,'geothermalFlux',0, ...
                'topMin',4,'topMax',4,'tunnelMin',4,'tunnelMax',4, ...
                'groundwaterInflowT',4,'groundwaterInflowGradient',0, ...
                'initial',struct('mode','profile','surfaceT',4,'gradient',0), ...
                'leakage',struct('jointConductance',zeros(1,8)), ...
                'thm',struct('enabled',true),'mechanics',struct('enable',false), ...
                'solver',struct('jacobian','analytic','reuseAcrossSteps',false));
            sourceUser=struct('outputFolder',sourceFolder,'maxYears',1,'model',model);
            source=run_TH_segment_year20_thm(sourceUser);
            sourceOut=source.lastYear;
            sourceArchive=load(fullfile(sourceFolder,'year_01_checkpoint.mat'));
            assert(size(sourceOut.cycles,1)==1&&sourceOut.timeDay==365, ...
                'THMCompletedTest:Fixture','The independent source must contain one complete year and cycle row.');

            % Direct public core recovery from a copied final checkpoint.
            coreFolder=fullfile(folder,'complete_core'); mkdir(coreFolder);
            copyFile=fullfile(coreFolder,'complete_checkpoint.mat');
            save(copyFile,'-struct','sourceArchive','-v7');
            coreModel=sourceOut.parameters;
            coreModel.initial=struct('mode','checkpoint','checkpointFile',copyFile,'resumeHistory',true);
            coreModel.outputPrefix=fullfile(coreFolder,'finalized');
            completeCore=TH_2D_segment_v16(coreModel);

            % A full year may have been checkpointed just before cycle/result
            % output. Preserve every original accepted record, remove only the
            % absent year diagnostic, and finalize the same full history.
            missingFolder=fullfile(folder,'missing_results');
            missingUser=TestTHMCompletedRecovery.cloneUncommitted(sourceFolder,sourceUser,missingFolder,false);
            monthArchive=load(fullfile(missingFolder,'year_01_checkpoint.mat'));
            monthArchive.cycles=zeros(0,5);
            save(fullfile(missingFolder,'year_01_checkpoint.mat'),'-struct','monthArchive','-v7');
            missingResults=run_TH_segment_year20_thm(missingUser);

            existingFolder=fullfile(folder,'existing_results');
            existingUser=TestTHMCompletedRecovery.cloneUncommitted(sourceFolder,sourceUser,existingFolder,true);
            existingResults=run_TH_segment_year20_thm(existingUser);

            alternateFolder=fullfile(folder,'unreadable_current');
            alternateUser=TestTHMCompletedRecovery.cloneUncommitted(sourceFolder,sourceUser,alternateFolder,true);
            currentCheckpointFile=fullfile(alternateFolder,'year_01_checkpoint.mat');
            alternateCheckpointFile=fullfile(alternateFolder,'year_01_FAILED_checkpoint.mat');
            copyfile(currentCheckpointFile,alternateCheckpointFile);
            fid=fopen(currentCheckpointFile,'wb');
            assert(fid>=0,'THMCompletedTest:Fixture','Cannot create the independent unreadable checkpoint fixture.');
            fwrite(fid,uint8('Independent deliberately unreadable checkpoint fixture.'));
            fclose(fid);
            alternateResults=run_TH_segment_year20_thm(alternateUser);
            r=struct('folder',folder,'sourceFolder',sourceFolder,'sourceUser',sourceUser, ...
                'sourceOut',sourceOut,'completeCore',completeCore, ...
                'missingResults',missingResults,'existingResults',existingResults, ...
                'alternateResults',alternateResults,'alternateCheckpointFile',alternateCheckpointFile);
            save(fullfile(folder,'completed_recovery_verification.mat'),'r','-v7');
            fprintf('Completed-history recovery verification folder: %s.\n',folder);
        end

        function user=cloneUncommitted(sourceFolder,sourceUser,folder,includeResults)
            mkdir(folder);
            copyfile(fullfile(sourceFolder,'year_01_checkpoint.mat'),fullfile(folder,'year_01_checkpoint.mat'));
            if includeResults
                copyfile(fullfile(sourceFolder,'year_01_results.mat'),fullfile(folder,'year_01_results.mat'));
            end
            d=load(fullfile(sourceFolder,'thm_progress.mat'),'progress'); progress=d.progress;
            % Simulate interruption before the driver's complete-year commit.
            progress.completedYears=0; progress.lastYear=[]; progress.coreCheckpoints={};
            progress.installedMechanicalReference=[];
            progress.annual=zeros(0,size(progress.annual,2));
            progress.monthly=zeros(0,size(progress.monthly,2));
            progress.byJoint=zeros(0,size(progress.byJoint,2));
            progress.daily=zeros(0,size(progress.daily,2));
            progress.config.outputFolder=folder;
            save(fullfile(folder,'thm_progress.mat'),'progress','-v7');
            user=sourceUser; user.outputFolder=folder; user.resume=true;
        end

        function corruptedArchive(r,kind)
            folder=tempname(r.folder);
            user=TestTHMCompletedRecovery.cloneUncommitted(r.sourceFolder,r.sourceUser,folder,true);
            file=fullfile(folder,'year_01_results.mat'); d=load(file,'out'); out=d.out;
            switch kind
                case 'physical'
                    out.parameters.Ksat=2*out.parameters.Ksat;
                case 'mesh'
                    out.mesh.pAll(1,1)=out.mesh.pAll(1,1)+1e-3;
                case 'version'
                    cpFile=fullfile(folder,'year_01_checkpoint.mat');
                    archive=load(cpFile);
                    archive.checkpoint.version='TH_v15_incompatible_fixture';
                    save(cpFile,'-struct','archive','-v7');
            end
            save(file,'out','-v7');
            run_TH_segment_year20_thm(user);
        end
    end
end
