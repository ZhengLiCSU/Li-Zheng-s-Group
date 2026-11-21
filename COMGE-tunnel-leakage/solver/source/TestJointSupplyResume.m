classdef TestJointSupplyResume < matlab.unittest.TestCase
    %TESTJOINTSUPPLYRESUME V14 monthly thermal/hydraulic history continuity.
    % Tests a genuine first-month checkpoint after a paused first year.
    % Each parameter uses a fresh folder and a warm periodic water table.
    % Existing full-history mechanism-study results are never overwritten.
    properties (TestParameter)
        freezeEnabled=struct('enabled',true,'disabled',false);
    end

    methods (Test)
        function monthlyResumePreservesHydraulicAndFreezeHistory(testCase,freezeEnabled)
            report=TestJointSupplyResume.exerciseResume(freezeEnabled);
            testCase.verifyEqual(report.pausedYears,1);
            testCase.verifyEqual(report.completedYears,2);
            testCase.verifyEqual(report.schemaVersion,8);
            testCase.verifyEqual(report.version,'TH_v14_20260930_joint_supply');
            testCase.verifyTrue(report.partialYearUsed);
            testCase.verifyTrue(report.phaseSnapshotsComplete);
            testCase.verifyEqual(report.finalTimeDay,730,'AbsTol',1e-7);
            testCase.verifyEqual(report.sumStepDays,365,'AbsTol',1e-6);
            testCase.verifyTrue(report.increasingDays);
            testCase.verifyEqual(report.mergedPrefixHistory,report.firstMonthHistory,'AbsTol',0);
            testCase.verifyEqual(report.mergedPrefixFreezeHistory,report.firstMonthFreezeHistory,'AbsTol',0);
            testCase.verifyEqual(report.mergedPrefixHydraulicHistory,report.firstMonthHydraulicHistory,'AbsTol',0);
            testCase.verifyEqual(report.partialPrefixHydraulicHistory,report.firstMonthHydraulicHistory,'AbsTol',0);
            testCase.verifyEqual(report.historyRows,report.freezeHistoryRows);
            testCase.verifyEqual(report.historyRows,report.hydraulicHistoryRows);
            testCase.verifyEqual(report.historyRows,report.acceptedSteps);
            testCase.verifyEqual(report.hydraulicHistoryColumns,27);
            testCase.verifyTrue(report.hydraulicRecordsValid);
            testCase.verifyEqual(report.actualDepthHistory,report.expectedDepthHistory,'AbsTol',1e-12);
            testCase.verifyGreaterThan(max(report.actualDepthHistory)-min(report.actualDepthHistory),.9);
            testCase.verifyEqual(report.dailyBins,365);
            testCase.verifyEqual(report.monthlyBins,12);
            testCase.verifyEqual(report.bodyVolume,0,'AbsTol',0);
            testCase.verifyEqual(report.totalVolume,report.jointVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.totalVolume,report.segmentVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.totalVolume,report.monthlyVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.totalVolume,report.dailyVolume,'AbsTol',1e-9);
            testCase.verifyEqual(report.sealedJointUnrestrictedVolume,zeros(1,7),'AbsTol',0);
            testCase.verifyEqual(report.freezeEnabled,freezeEnabled);
            testCase.verifyTrue(report.temperatureRecordsValid);
            testCase.verifyTrue(report.factorRecordsValid);
            testCase.verifyTrue(report.year20ResponseIsAbsent);
            testCase.verifyLessThan(report.maximumWaterDefect_m2,1e-6);
            testCase.verifyLessThan(report.maximumEnergyDefect_J_per_m,100);
        end
    end

    methods (Static,Access=private)
        function report=exerciseResume(freezeEnabled)
            root=fileparts(mfilename('fullpath')); addpath(root);
            suffix='disabled';
            if freezeEnabled, suffix='enabled'; end
            parent=fullfile(root,'joint_supply_study_20260930',['resume_check_' suffix]);
            if ~isfolder(parent), mkdir(parent); end
            folder=tempname(parent); mkdir(folder);
            base=joint_only_th_model();
            day=(0:365)'; depth=5-.5*cos(2*pi*(day-135)/365); depth(end)=depth(1);
            model=struct('dtMax',4*86400, ...
                'topMin',base.topMin+10,'topMax',base.topMax+10, ...
                'tunnelMin',base.tunnelMin+10,'tunnelMax',base.tunnelMax+10, ...
                'groundwaterInflowT',base.groundwaterInflowT+10, ...
                'initial',struct('surfaceT',base.initial.surfaceT+10), ...
                'forcing',struct('day',day,'waterTableDepth_m',depth), ...
                'leakage',struct('jointFreeze',struct('enabled',freezeEnabled)));
            user=struct('model',model,'maxYears',2,'consecutiveYears',21, ...
                'productionYears',1,'outputFolder',folder);
            marker=fullfile(folder,'PAUSE_AFTER_YEAR');
            fid=fopen(marker,'w'); assert(fid>=0,'Cannot create the staged pause marker.'); fclose(fid);
            paused=run_TH_segment_year20_regime(user);
            assert(paused.completedYears==1,'The staged test did not pause after year one.');
            delete(marker);

            first=paused.lastYear;
            opt=first.parameters; opt.mesh=first.mesh;
            opt.initial=struct('mode','fields','T',first.T,'h',first.h,'startDay',365);
            opt.durationDays=365/12; opt.outputEveryDays=365/12;
            opt.saveOutput=true; opt.plot=false; opt.verbose=false; opt.mechanics.enable=false;
            opt.outputPrefix=fullfile(folder,'active_year');
            month=TH_2D_segment_v14(opt);
            assert(abs(month.timeDay-(365+365/12))<1e-7,'The first-month checkpoint is incomplete.');
            user.resume=true;
            resumed=run_TH_segment_year20_regime(user);
            last=resumed.lastYear; leak=last.leakage;
            freeze=leak.jointFreeze; hydraulic=leak.hydraulics;
            prefixRows=size(month.history,1);
            prefix=load(fullfile(folder,'partial_year_prefix.mat'),'partialYearBase');
            if freezeEnabled
                temperaturesValid=all(isfinite(freeze.minTemperature_C(:))) ...
                    &&all(isfinite(freeze.maxTemperature_C(:))) ...
                    &&all(freeze.minTemperature_C(:)<=freeze.maxTemperature_C(:));
                factorsValid=all(isfinite(freeze.factor(:)))&&all(freeze.factor(:)>=0);
            else
                temperaturesValid=all(isnan(freeze.minTemperature_C(:))) ...
                    &&all(isnan(freeze.maxTemperature_C(:)));
                factorsValid=all(freeze.factor(:)==1)&&all(freeze.frozenLengthFraction(:)==0);
            end
            hydraulicValid=all(isfinite(hydraulic.history(:))) ...
                &&all(hydraulic.unrestrictedRate_m2_per_s(:)>=0) ...
                &&all(hydraulic.negativeTransEdges>=0) ...
                &&all(hydraulic.negativeTransEdges==fix(hydraulic.negativeTransEdges));
            phase=mod(last.history(:,1),365);
            expectedDepth=interp1(day,depth,phase,'linear');
            dt=last.history(:,2)*86400;
            unrestrictedVolume=dt'*hydraulic.unrestrictedRate_m2_per_s;
            sealed=last.parameters.leakage.jointConductance==0;
            report=struct('outputFolder',folder,'pausedYears',paused.completedYears, ...
                'completedYears',resumed.completedYears,'schemaVersion',resumed.schemaVersion, ...
                'version',last.version,'partialYearUsed',last.partialYearResume.used, ...
                'phaseSnapshotsComplete',last.partialYearResume.phaseSnapshotsComplete, ...
                'finalTimeDay',last.timeDay,'sumStepDays',sum(last.history(:,2)), ...
                'increasingDays',all(diff(last.history(:,1))>0), ...
                'mergedPrefixHistory',last.history(1:prefixRows,:),'firstMonthHistory',month.history, ...
                'mergedPrefixFreezeHistory',freeze.history(1:prefixRows,:), ...
                'firstMonthFreezeHistory',month.leakage.jointFreeze.history, ...
                'mergedPrefixHydraulicHistory',hydraulic.history(1:prefixRows,:), ...
                'firstMonthHydraulicHistory',month.leakage.hydraulics.history, ...
                'partialPrefixHydraulicHistory',prefix.partialYearBase.jointHydraulicHist, ...
                'historyRows',size(last.history,1),'freezeHistoryRows',size(freeze.history,1), ...
                'hydraulicHistoryRows',size(hydraulic.history,1), ...
                'hydraulicHistoryColumns',size(hydraulic.history,2), ...
                'acceptedSteps',last.solverWork.acceptedSteps,'hydraulicRecordsValid',hydraulicValid, ...
                'actualDepthHistory',hydraulic.waterTableDepth_m,'expectedDepthHistory',expectedDepth, ...
                'dailyBins',numel(leak.dailyTotalVolume_m3_per_m), ...
                'monthlyBins',numel(leak.monthlyTotalVolume_m3_per_m), ...
                'bodyVolume',leak.bodyVolume_m3_per_m,'totalVolume',leak.totalVolume_m3_per_m, ...
                'jointVolume',sum(leak.volumeByJoint_m3_per_m), ...
                'segmentVolume',sum(leak.volumeBySegment_m3_per_m), ...
                'monthlyVolume',sum(leak.monthlyTotalVolume_m3_per_m), ...
                'dailyVolume',sum(leak.dailyTotalVolume_m3_per_m), ...
                'sealedJointUnrestrictedVolume',unrestrictedVolume(sealed), ...
                'freezeEnabled',freeze.enabled,'temperatureRecordsValid',temperaturesValid, ...
                'factorRecordsValid',factorsValid, ...
                'year20ResponseIsAbsent',isnan(resumed.year20AnnualLeakage_L_per_m), ...
                'maximumWaterDefect_m2',max(abs(last.history(:,8))), ...
                'maximumEnergyDefect_J_per_m',max(abs(last.history(:,9))));
            save(fullfile(folder,'resume_verification_report.mat'),'report','-v7');
            fprintf('V14 resume verification: enabled=%d, folder=%s.\n',freezeEnabled,folder);
        end
    end
end
