classdef TestJointFreezeResume < matlab.unittest.TestCase
    %TESTJOINTFREEZERESUME Preserve complete year-2 leakage/icing on restart.
    % Each parameter creates a fresh subfolder. Existing full-20-year case
    % files and earlier test records are neither read nor overwritten.
    % Run results=runtests('TestJointFreezeResume.m'); assertSuccess(results).

    properties (TestParameter)
        freezeEnabled=struct('enabled',true,'disabled',false);
    end

    methods (Test)
        function warmMonthlyResumePreservesYearHistory(testCase,freezeEnabled)
            % Arrange and act: the utility contains the staged numerical run.
            report=TestJointFreezeResume.exerciseResume(freezeEnabled);

            % Assert: saved first-month records survive the annual merge.
            testCase.verifyEqual(report.pausedYears,1);
            testCase.verifyEqual(report.completedYears,2);
            testCase.verifyEqual(report.schemaVersion,7);
            testCase.verifyTrue(report.partialYearUsed);
            testCase.verifyEqual(report.finalTimeDay,730,'AbsTol',1e-7);
            testCase.verifyEqual(report.sumStepDays,365,'AbsTol',1e-6);
            testCase.verifyTrue(report.increasingDays);
            testCase.verifyEqual(report.mergedPrefixHistory, ...
                report.firstMonthHistory,'AbsTol',0);
            testCase.verifyEqual(report.mergedPrefixJointFreezeHistory, ...
                report.firstMonthJointFreezeHistory,'AbsTol',0);
            testCase.verifyEqual(report.historyRows,report.freezeHistoryRows);
            testCase.verifyEqual(report.historyRows,report.acceptedSteps);
            testCase.verifyEqual(report.dailyBins,365);
            testCase.verifyEqual(report.monthlyBins,12);
            testCase.verifyEqual(report.bodyVolume,0,'AbsTol',0);
            testCase.verifyEqual(report.totalVolume,report.jointVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.totalVolume,report.segmentVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.totalVolume,report.monthlyVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.totalVolume,report.dailyVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.freezeEnabled,freezeEnabled);
            testCase.verifyTrue(report.temperatureRecordsValid);
            testCase.verifyTrue(report.factorRecordsValid);
            testCase.verifyLessThan(report.maximumWaterDefect_m2,1e-6);
            testCase.verifyLessThan(report.maximumEnergyDefect_J_per_m,100);
        end
    end

    methods (Static,Access=private)
        function report=exerciseResume(freezeEnabled)
            root=fileparts(mfilename('fullpath'));
            addpath(root);
            if freezeEnabled
                parent=fullfile(root,'joint_freeze_study_20260930','resume_check_enabled');
            else
                parent=fullfile(root,'joint_freeze_study_20260930','resume_check_disabled');
            end
            if ~isfolder(parent), mkdir(parent); end
            folder=tempname(parent);
            mkdir(folder);
            base=joint_only_th_model();
            model=struct('dtMax',4*86400, ...
                'topMin',base.topMin+10,'topMax',base.topMax+10, ...
                'tunnelMin',base.tunnelMin+10,'tunnelMax',base.tunnelMax+10, ...
                'groundwaterInflowT',base.groundwaterInflowT+10, ...
                'initial',struct('surfaceT',base.initial.surfaceT+10), ...
                'leakage',struct('jointFreeze',struct('enabled',freezeEnabled)));
            user=struct('model',model,'maxYears',2,'consecutiveYears',21, ...
                'productionYears',1,'outputFolder',folder);
            marker=fullfile(folder,'PAUSE_AFTER_YEAR');
            fid=fopen(marker,'w');
            assert(fid>0,'Cannot create the test pause marker.');
            fclose(fid);
            paused=run_TH_segment_year20_jointfreeze(user);
            assert(paused.completedYears==1, ...
                'The staged test must pause after its first complete year.');
            delete(marker);

            % Create a genuine core month checkpoint in year 2, preserving
            % all model parameters while replacing only time/initial state.
            first=paused.lastYear;
            opt=first.parameters;
            opt.mesh=first.mesh;
            opt.initial=struct('mode','fields','T',first.T,'h',first.h,'startDay',365);
            opt.durationDays=365/12; opt.outputEveryDays=365/12;
            opt.saveOutput=true; opt.plot=false; opt.verbose=false;
            opt.mechanics.enable=false;
            opt.outputPrefix=fullfile(folder,'active_year');
            month=TH_2D_segment_v12(opt);
            assert(abs(month.timeDay-(365+365/12))<1e-7, ...
                'The first month checkpoint did not reach its report boundary.');
            user.resume=true;
            resumed=run_TH_segment_year20_jointfreeze(user);
            last=resumed.lastYear;
            leak=last.leakage;
            freeze=leak.jointFreeze;
            prefixRows=size(month.history,1);
            if freezeEnabled
                temperaturesValid=all(isfinite(freeze.minTemperature_C(:))) && ...
                    all(isfinite(freeze.maxTemperature_C(:))) && ...
                    all(freeze.minTemperature_C(:)<=freeze.maxTemperature_C(:));
                factorsValid=all(isfinite(freeze.factor(:))) && all(freeze.factor(:)>=0);
            else
                temperaturesValid=all(isnan(freeze.minTemperature_C(:))) && ...
                    all(isnan(freeze.maxTemperature_C(:)));
                factorsValid=all(freeze.factor(:)==1) && ...
                    all(freeze.frozenLengthFraction(:)==0);
            end
            report=struct('outputFolder',folder, ...
                'pausedYears',paused.completedYears,'completedYears',resumed.completedYears, ...
                'schemaVersion',resumed.schemaVersion, ...
                'partialYearUsed',last.partialYearResume.used,'finalTimeDay',last.timeDay, ...
                'sumStepDays',sum(last.history(:,2)), ...
                'increasingDays',all(diff(last.history(:,1))>0), ...
                'mergedPrefixHistory',last.history(1:prefixRows,:), ...
                'firstMonthHistory',month.history, ...
                'mergedPrefixJointFreezeHistory',freeze.history(1:prefixRows,:), ...
                'firstMonthJointFreezeHistory',month.leakage.jointFreeze.history, ...
                'historyRows',size(last.history,1),'freezeHistoryRows',size(freeze.history,1), ...
                'acceptedSteps',last.solverWork.acceptedSteps, ...
                'dailyBins',numel(leak.dailyTotalVolume_m3_per_m), ...
                'monthlyBins',numel(leak.monthlyTotalVolume_m3_per_m), ...
                'bodyVolume',leak.bodyVolume_m3_per_m, ...
                'totalVolume',leak.totalVolume_m3_per_m, ...
                'jointVolume',sum(leak.volumeByJoint_m3_per_m), ...
                'segmentVolume',sum(leak.volumeBySegment_m3_per_m), ...
                'monthlyVolume',sum(leak.monthlyTotalVolume_m3_per_m), ...
                'dailyVolume',sum(leak.dailyTotalVolume_m3_per_m), ...
                'freezeEnabled',freeze.enabled, ...
                'temperatureRecordsValid',temperaturesValid,'factorRecordsValid',factorsValid, ...
                'maximumWaterDefect_m2',max(abs(last.history(:,8))), ...
                'maximumEnergyDefect_J_per_m',max(abs(last.history(:,9))));
            save(fullfile(folder,'resume_verification_report.mat'),'report','-v7');
            fprintf('Resume verification: enabled=%d, folder=%s.\n',freezeEnabled,folder);
        end
    end
end
