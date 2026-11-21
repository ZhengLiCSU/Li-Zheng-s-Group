classdef TestTHMPartialResume < matlab.unittest.TestCase
    %TESTTHMPARTIALRESUME A month-style checkpoint resumes an original run.
    % Small warm, saturated leaking annulus; this is a restart integration
    % test, not the research mesh or a freeze-thaw field calibration. The
    % 0.05-day checkpoint is reconstructed from the reference 0.10-day
    % accepted prefix and its actual snapshot. A separate matching prefix
    % run supplies the true prefix local-ledger maxima; future full-run
    % maxima are never inserted into the synthetic earlier checkpoint.
    % Three public core calls are made once in class setup. No core/driver
    % edits, timer interruption or unregistered partial-year concatenation.

    properties (Access=private)
        Evidence
    end

    methods (TestClassSetup)
        function prepareActualRestartEvidence(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.Evidence=TestTHMPartialResume.exercise(root);
        end
    end

    methods (Test)
        function entireStateAndOriginalRunStartArePreserved(testCase)
            r=testCase.Evidence;
            testCase.verifyEqual(r.resumed.finalX,r.reference.finalX,'AbsTol',2e-9);
            testCase.verifyEqual(r.resumed.initialX,r.reference.initialX,'AbsTol',0);
            testCase.verifyEqual(r.resumed.timeDay,.1,'AbsTol',1e-12);
            testCase.verifyEqual(sum(r.resumed.history(:,2)),.1,'AbsTol',1e-12);
            testCase.verifyEqual(r.partialCheckpoint.time/86400,.05,'AbsTol',1e-12);
            testCase.verifyEqual(r.partialCheckpoint.runInventory.startTime,0,'AbsTol',0);
            testCase.verifyEqual(r.partialCheckpoint.runInventory.endTime/86400,.1,'AbsTol',1e-12);
            testCase.verifyEqual(r.finalCheckpoint.X,r.resumed.finalX,'AbsTol',0);
        end

        function everyAcceptedDateAndPrefixIsUnique(testCase)
            r=testCase.Evidence;
            testCase.verifyEqual(r.resumed.history(:,1:2),r.reference.history(:,1:2),'AbsTol',0);
            testCase.verifyEqual(r.resumed.history(1:r.prefixCount,:),r.reference.history(1:r.prefixCount,:),'AbsTol',0);
            testCase.verifyEqual(r.resumed.history(:,1),[.025;.05;.075;.1],'AbsTol',1e-12);
            testCase.verifyEqual(numel(unique(r.resumed.history(:,1))),size(r.resumed.history,1));
            testCase.verifyEqual(r.resumed.solverWork.acceptedSteps,r.reference.solverWork.acceptedSteps);
            testCase.verifyEqual(r.resumed.solverWork.rejectedSteps,r.reference.solverWork.rejectedSteps);
            testCase.verifyEqual(r.finalCheckpointArchive.hist,r.resumed.history,'AbsTol',0);
        end

        function conservativeJointDailyAndCalendarMonthBudgetsMatch(testCase)
            r=testCase.Evidence;
            testCase.verifyGreaterThan(r.reference.leakage.totalVolume_m3_per_m,0);
            testCase.verifyEqual(r.resumed.leakage.jointRate_m2_per_s,r.reference.leakage.jointRate_m2_per_s,'AbsTol',1e-18);
            testCase.verifyEqual(r.resumed.leakage.segmentRate_m2_per_s,r.reference.leakage.segmentRate_m2_per_s,'AbsTol',1e-18);
            testCase.verifyEqual(r.resumed.leakage.volumeByJoint_m3_per_m,r.reference.leakage.volumeByJoint_m3_per_m,'AbsTol',1e-12);
            testCase.verifyEqual(r.resumed.leakage.totalVolume_m3_per_m,r.reference.leakage.totalVolume_m3_per_m,'AbsTol',1e-12);
            testCase.verifyEqual(r.resumed.leakage.dailyBinEndDay,r.reference.leakage.dailyBinEndDay,'AbsTol',0);
            testCase.verifyEqual(r.resumed.leakage.dailyVolumeByJoint_m3_per_m,r.reference.leakage.dailyVolumeByJoint_m3_per_m,'AbsTol',1e-12);
            testCase.verifyEqual(r.resumedPartialJanuaryJointVolume,r.referencePartialJanuaryJointVolume,'AbsTol',1e-12);
            testCase.verifyEqual(r.resumedPartialJanuaryJointVolume,r.resumed.leakage.volumeByJoint_m3_per_m,'AbsTol',1e-12);
            testCase.verifyEqual(r.resumed.leakage.bodyVolume_m3_per_m,0,'AbsTol',0);
            testCase.verifyEqual(r.waterAnnualSumDifference_m2,0,'AbsTol',1e-10);
            testCase.verifyEqual(r.energyAnnualSumDifference_J_per_m,0,'AbsTol',.01);
        end

        function mechanicalAndHydraulicHistoriesKeepTheCompletePrefix(testCase)
            r=testCase.Evidence;
            testCase.verifyEqual(r.resumed.mechanics.history(1:r.prefixCount,:),r.reference.mechanics.history(1:r.prefixCount,:),'AbsTol',0);
            testCase.verifyEqual(r.resumed.mechanics.jointHistory(1:r.prefixCount,:),r.reference.mechanics.jointHistory(1:r.prefixCount,:),'AbsTol',0);
            testCase.verifyEqual(r.resumed.mechanics.strainHistory(1:r.prefixCount,:),r.reference.mechanics.strainHistory(1:r.prefixCount,:),'AbsTol',0);
            testCase.verifyEqual(r.resumed.leakage.hydraulics.history(1:r.prefixCount,:),r.reference.leakage.hydraulics.history(1:r.prefixCount,:),'AbsTol',0);
            testCase.verifyEqual(r.finalCheckpointArchive.mechanicalHistory,r.resumed.mechanics.history,'AbsTol',0);
            testCase.verifyEqual(r.finalCheckpointArchive.jointMechanicalHistory,r.resumed.mechanics.jointHistory,'AbsTol',0);
            testCase.verifyEqual(r.finalCheckpointArchive.strainHistory,r.resumed.mechanics.strainHistory,'AbsTol',0);
            testCase.verifyEqual(r.finalCheckpointArchive.jointHydraulicHist,r.resumed.leakage.hydraulics.history,'AbsTol',0);
        end

        function actualInstalledReferenceAndSnapshotsSurviveRestart(testCase)
            r=testCase.Evidence;
            testCase.verifyEqual(r.resumed.mechanics.referenceU,r.reference.mechanics.referenceU,'AbsTol',0);
            testCase.verifyEqual(r.resumed.mechanics.referenceTemperature,r.reference.mechanics.referenceTemperature,'AbsTol',0);
            testCase.verifyEqual(r.finalCheckpoint.mechanicalReference,r.partialCheckpoint.mechanicalReference);
            testCase.verifyEqual([r.resumed.snapshots.timeDay],[r.reference.snapshots.timeDay],'AbsTol',0);
            testCase.verifyEqual(r.resumed.snapshots(1).U,r.reference.snapshots(1).U,'AbsTol',0);
            testCase.verifyEqual(r.resumed.snapshots(end).U,r.reference.snapshots(end).U,'AbsTol',2e-9);
            testCase.verifyEqual(r.resumed.version,'THM_v16_20260930_analytic_p1');
            testCase.verifyTrue(r.resumed.parameters.initial.resumeHistory);
        end

        function prefixLocalLedgerInventoryUsesNoFutureInformation(testCase)
            r=testCase.Evidence;
            testCase.verifyEqual(r.partialCheckpoint.runInventory.maxLocalWaterDefect, ...
                r.prefix.numerics.maximumAcceptedLocalWaterDefect_m2,'AbsTol',0);
            testCase.verifyEqual(r.partialCheckpoint.runInventory.maxLocalEnergyDefect, ...
                r.prefix.numerics.maximumAcceptedLocalEnergyDefect_Jm,'AbsTol',0);
            testCase.verifyEqual(r.resumed.numerics.maximumAcceptedLocalWaterDefect_m2, ...
                r.reference.numerics.maximumAcceptedLocalWaterDefect_m2,'AbsTol',1e-11);
            testCase.verifyEqual(r.resumed.numerics.maximumAcceptedLocalEnergyDefect_Jm, ...
                r.reference.numerics.maximumAcceptedLocalEnergyDefect_Jm,'AbsTol',.01);
        end
    end

    methods (Static,Access=private)
        function r=exercise(root)
            parent=fullfile(root,'thm_partial_resume_checks');
            if ~isfolder(parent),mkdir(parent);end
            folder=tempname(parent);mkdir(folder);
            model=joint_contact_thm_model(struct('mesh',make_thm_verification_mesh(16), ...
                'durationDays',.1,'outputEveryDays',.05,'dtInitial',.025*86400, ...
                'dtMax',.025*86400,'dtMin',.025*86400/128, ...
                'waterTableDepth',-2,'topSeepage',false,'geothermalFlux',0, ...
                'topMin',4,'topMax',4,'tunnelMin',4,'tunnelMax',4, ...
                'groundwaterInflowT',4,'groundwaterInflowGradient',0, ...
                'initial',struct('mode','profile','surfaceT',4,'gradient',0,'startDay',0), ...
                'mechanics',struct('enable',false),'thm',struct('enabled',true), ...
                'solver',struct('tol',1e-10,'maxIter',60,'predictor',false, ...
                    'reuseAcrossSteps',false,'jacobian','analytic'), ...
                'plot',false,'verbose',false,'saveOutput',true, ...
                'outputPrefix',fullfile(folder,'continuous')));
            reference=TH_2D_segment_v16(model);
            source=load([model.outputPrefix '_checkpoint.mat']);
            index=find(reference.history(:,1)<=.05+1e-12);
            snapshotIndex=find(abs([reference.snapshots.timeDay]-.05)<1e-12);
            assert(numel(index)==2&&numel(snapshotIndex)==1,'THMPartialTest:Fixture', ...
                'The reference must have two accepted prefix steps and a real 0.05-day snapshot.');
            snapshot=reference.snapshots(snapshotIndex);
            prefixModel=model;prefixModel.durationDays=.05;prefixModel.outputEveryDays=.05;
            prefixModel.outputPrefix=fullfile(folder,'prefix_inventory');
            prefix=TH_2D_segment_v16(prefixModel);
            partialX=[snapshot.T;snapshot.h;snapshot.U];
            assert(isequal(partialX,prefix.finalX)&&isequal(prefix.initialX,reference.initialX), ...
                'THMPartialTest:Fixture','Independent prefix run must match the original snapshot exactly.');
            assert(isequal(prefix.history,reference.history(index,:)), ...
                'THMPartialTest:Fixture','The independent prefix run must follow identical accepted steps.');

            % Only original accepted prefix records enter the recovered MAT.
            source.checkpoint.X=partialX;source.checkpoint.time=.05*86400;
            source.checkpoint.rejectedSteps=prefix.solverWork.rejectedSteps;
            source.checkpoint.runInventory.startX=reference.initialX;
            source.checkpoint.runInventory.startTime=0;
            source.checkpoint.runInventory.endTime=.1*86400;
            source.checkpoint.runInventory.nextReport=.1*86400;
            source.checkpoint.runInventory.nextYear=365*86400;
            source.checkpoint.runInventory.yearRef=reference.initialX;
            source.checkpoint.runInventory.nextDt=.025*86400;
            source.checkpoint.runInventory.priorStart=[];
            source.checkpoint.runInventory.priorAcceptedDt=NaN;
            source.checkpoint.runInventory.maxLocalWaterDefect=prefix.numerics.maximumAcceptedLocalWaterDefect_m2;
            source.checkpoint.runInventory.maxLocalEnergyDefect=prefix.numerics.maximumAcceptedLocalEnergyDefect_Jm;
            source.hist=reference.history(index,:);source.cycles=zeros(0,5);
            source.snaps=reference.snapshots(snapshotIndex);
            source.leakHistSegment=reference.leakage.segmentRate_m2_per_s(index,:);
            source.leakHistJoint=reference.leakage.jointRate_m2_per_s(index,:);
            source.leakHistBody=reference.leakage.bodyRate_m2_per_s(index,:);
            source.jointFreezeHist=source.jointFreezeHist(index,:);
            source.jointHydraulicHist=reference.leakage.hydraulics.history(index,:);
            source.mechanicalHistory=reference.mechanics.history(index,:);
            source.jointMechanicalHistory=reference.mechanics.jointHistory(index,:);
            source.strainHistory=reference.mechanics.strainHistory(index,:);
            partialFile=fullfile(folder,'month_style_checkpoint.mat');
            save(partialFile,'-struct','source','-v7');
            resumeModel=model;resumeModel.outputPrefix=fullfile(folder,'resumed');
            resumeModel.initial=struct('mode','checkpoint','checkpointFile',partialFile,'resumeHistory',true);
            resumed=TH_2D_segment_v16(resumeModel);
            finalArchive=load([resumeModel.outputPrefix '_checkpoint.mat']);
            % This sub-day verification occupies only a PARTIAL January;
            % it cannot claim twelve complete month bins or a full year.
            JanuaryReference=TestTHMPartialResume.partialJanuaryJointVolume(reference);
            JanuaryResumed=TestTHMPartialResume.partialJanuaryJointVolume(resumed);
            r=struct('outputFolder',folder,'reference',reference,'prefix',prefix,'resumed',resumed, ...
                'prefixCount',numel(index),'partialCheckpoint',source.checkpoint, ...
                'finalCheckpoint',finalArchive.checkpoint,'finalCheckpointArchive',finalArchive, ...
                'referencePartialJanuaryJointVolume',JanuaryReference, ...
                'resumedPartialJanuaryJointVolume',JanuaryResumed, ...
                'waterAnnualSumDifference_m2',sum(resumed.history(:,8))-sum(reference.history(:,8)), ...
                'energyAnnualSumDifference_J_per_m',sum(resumed.history(:,9))-sum(reference.history(:,9)));
            save(fullfile(folder,'partial_resume_verification.mat'),'r','-v7');
            fprintf('THM partial restart verification folder: %s.\n',folder);
        end

        function volume=partialJanuaryJointVolume(out)
            finish=out.history(:,1);start=finish-out.history(:,2);
            overlap=max(0,min(finish,31)-max(start,0));
            volume=sum(out.leakage.jointRate_m2_per_s.*(overlap*86400),1);
        end
    end
end
