function study=verify_joint_freeze_year20(folder)
%VERIFY_JOINT_FREEZE_YEAR20 Last-year time/path integration sensitivity.
% This numerical research artifact reruns ONLY year 20 from the identical
% end-of-year-19 coldIce state. It is not full-20-year time convergence,
% spatial mesh convergence, joint-law calibration, or a unit test suite.
% The existing full-history dt4/nq8 result is compared with dt1/nq8 and
% dt0.5/nq8; path samples nq4 and nq16 are compared with nq8 at dt1.
% Each new case writes MAT, accepted-step CSV, daily-volume CSV and a
% summary CSV. No probability/physics claim is based on error thresholds.
if nargin<1 || isempty(folder)
    root=fileparts(mfilename('fullpath'));
    folder=fullfile(root,'joint_freeze_study_20260930');
end
root=fileparts(mfilename('fullpath')); addpath(root);
sourceFile=fullfile(folder,'coldIce','TH_segment_spinup_results.mat');
assert(isfile(sourceFile), ...
    'Complete coldIce full-history result is required before this study.');
loaded=load(sourceFile,'result'); r=loaded.result;
assert(r.completedYears==20 && ...
    strcmp(r.lastYear.version,'TH_v12_20260930_joint_freeze'), ...
    'The input must be a complete v12 twenty-year coldIce result.');
base=r.lastYear;
assert(base.parameters.leakage.jointFreeze.enabled && ...
    base.parameters.leakage.jointFreeze.iceEnabled, ...
    'The baseline must have reversible joint icing enabled.');
assert(base.parameters.leakage.jointFreeze.quadraturePoints==8 && ...
    abs(base.parameters.dtMax/86400-4)<1e-12, ...
    'This comparison expects the existing dt4/nq8 coldIce baseline.');
thermalNodes=numel(base.T);
assert(numel(base.initialX)==thermalNodes+numel(base.h), ...
    'The year-19 final state is missing or has an incompatible shape.');
assert(abs(base.timeDay-20*365)<1e-7 && ...
    abs(base.history(1,1)-base.history(1,2)-19*365)<1e-7, ...
    'Baseline history must represent the complete twentieth year.');
outputFolder=fullfile(folder,'year20_time_path_check');
if ~isfolder(outputFolder), mkdir(outputFolder); end

caseName={'existing_dt4_nq8';'dt1_nq8';'dtHalf_nq8';'dt1_nq4';'dt1_nq16'};
dtMaxDays=[4;1;.5;1;1];
pathPoints=[8;8;8;4;16];
outputs=cell(5,1); outputs{1}=base;
elapsed_s=nan(5,1);
studyDefinition=struct('sourceFile',sourceFile, ...
    'initialStateSource','coldIce.lastYear.initialX: end of year 19', ...
    'startDay',19*365,'durationDays',365, ...
    'timeReferenceCase','dtHalf_nq8', ...
    'pathReferenceCase','dt1_nq8', ...
    'limitation',['Only final-year integration and path quadrature are examined. ' ...
        'The prior nineteen-year state and spatial mesh are fixed. This is not ' ...
        'full-history time convergence or spatial mesh convergence.'], ...
    'baselineLastYearRuntimeStatus',['Not separately recorded; NaN in the table. ' ...
        'Complete twenty-year wall time is not used as final-year runtime.']);
for index=2:5
    opt=base.parameters;
    opt.mesh=base.mesh;
    opt.initial=struct('mode','fields', ...
        'T',base.initialX(1:thermalNodes), ...
        'h',base.initialX(thermalNodes+1:end),'startDay',19*365);
    opt.durationDays=365; opt.outputEveryDays=365/12;
    opt.dtInitial=.1*86400;
    opt.dtMax=dtMaxDays(index)*86400;
    opt.dtMin=86400/16384;
    opt.leakage.jointFreeze.quadraturePoints=pathPoints(index);
    opt.mechanics.enable=false;
    opt.saveOutput=true; opt.plot=false; opt.verbose=false;
    opt.outputPrefix=fullfile(outputFolder,caseName{index});
    fprintf('Year-20 integration study: %s (same year-19 state).\n',caseName{index});
    started=tic;
    out=TH_2D_segment_v12(opt);
    elapsed_s(index)=toc(started);
    assertCompleteYear(out);
    outputs{index}=out;
    metadata=struct('caseName',caseName{index},'elapsed_s',elapsed_s(index), ...
        'dtMaxDays',dtMaxDays(index),'pathPoints',pathPoints(index), ...
        'studyDefinition',studyDefinition); %#ok<NASGU>
    save([opt.outputPrefix '_results.mat'],'metadata','-append');
    exportHistory(out,outputFolder,caseName{index});
end

Q20_L_per_m=zeros(5,1); acceptedSteps=zeros(5,1); rejectedSteps=zeros(5,1);
sumAbsWaterDefect_m2=zeros(5,1); sumAbsEnergyDefect_J_per_m=zeros(5,1);
maxWaterDefect_m2=zeros(5,1); maxEnergyDefect_J_per_m=zeros(5,1);
sumWaterDefect_m2=zeros(5,1); sumEnergyDefect_J_per_m=zeros(5,1);
minimumStepDays=zeros(5,1); maximumResidual=zeros(5,1);
dailyVolume=cell(5,1); monthlyVolume=cell(5,1);
for index=1:5
    out=outputs{index}; h=out.history;
    Q20_L_per_m(index)=1000*out.leakage.totalVolume_m3_per_m;
    acceptedSteps(index)=out.solverWork.acceptedSteps;
    rejectedSteps(index)=out.solverWork.rejectedSteps;
    sumAbsWaterDefect_m2(index)=sum(abs(h(:,8)));
    sumAbsEnergyDefect_J_per_m(index)=sum(abs(h(:,9)));
    maxWaterDefect_m2(index)=max(abs(h(:,8)));
    maxEnergyDefect_J_per_m(index)=max(abs(h(:,9)));
    sumWaterDefect_m2(index)=sum(h(:,8));
    sumEnergyDefect_J_per_m(index)=sum(h(:,9));
    minimumStepDays(index)=min(h(:,2));
    maximumResidual(index)=max(h(:,11));
    dailyVolume{index}=1000*out.leakage.dailyTotalVolume_m3_per_m(:);
    monthlyVolume{index}=1000*out.leakage.monthlyTotalVolume_m3_per_m(:);
end
relativeToDtHalf_nq8_percent=100*abs(Q20_L_per_m-Q20_L_per_m(3)) ...
    /max(Q20_L_per_m(3),realmin);
relativeToDt1_nq8_percent=100*abs(Q20_L_per_m-Q20_L_per_m(2)) ...
    /max(Q20_L_per_m(2),realmin);
referenceCase={'dtHalf_nq8';'dtHalf_nq8';'dtHalf_nq8';'dt1_nq8';'dt1_nq8'};
primaryRelativeDifference_percent=[relativeToDtHalf_nq8_percent(1:3); ...
    relativeToDt1_nq8_percent(4:5)];
dailyL1Difference_L_per_m=zeros(5,1);
maximumMonthlyDifference_L_per_m=zeros(5,1);
referenceIndex=[3;3;3;2;2];
for index=1:5
    reference=referenceIndex(index);
    dailyL1Difference_L_per_m(index)=sum(abs(dailyVolume{index}-dailyVolume{reference}));
    maximumMonthlyDifference_L_per_m(index)=max(abs( ...
        monthlyVolume{index}-monthlyVolume{reference}));
end
dailyL1RelativeDifference_percent=100*dailyL1Difference_L_per_m ...
    ./max(Q20_L_per_m(referenceIndex),realmin);
comparison=table(caseName,dtMaxDays,pathPoints,Q20_L_per_m, ...
    relativeToDtHalf_nq8_percent,relativeToDt1_nq8_percent,referenceCase, ...
    primaryRelativeDifference_percent,dailyL1Difference_L_per_m, ...
    dailyL1RelativeDifference_percent,maximumMonthlyDifference_L_per_m, ...
    acceptedSteps,rejectedSteps,elapsed_s,minimumStepDays,maximumResidual, ...
    sumAbsWaterDefect_m2,sumAbsEnergyDefect_J_per_m,maxWaterDefect_m2, ...
    maxEnergyDefect_J_per_m,sumWaterDefect_m2,sumEnergyDefect_J_per_m);
writetable(comparison,fullfile(outputFolder,'year20_time_path_comparison.csv'));
for index=2:5
    caseSummary=comparison(index,:); %#ok<NASGU>
    writetable(caseSummary,fullfile(outputFolder,[caseName{index} '_summary.csv']));
    save(fullfile(outputFolder,[caseName{index} '_results.mat']), ...
        'caseSummary','studyDefinition','-append');
end
study=struct('definition',studyDefinition,'comparison',comparison, ...
    'dailyVolume_L_per_m',{dailyVolume},'monthlyVolume_L_per_m',{monthlyVolume});
save(fullfile(outputFolder,'year20_time_path_study.mat'),'study','-v7');
disp(comparison);
fprintf('%s\n',studyDefinition.limitation);
end

function assertCompleteYear(out)
leak=out.leakage;
assert(abs(out.timeDay-20*365)<1e-7 && ...
    abs(sum(out.history(:,2))-365)<1e-6 && all(isfinite(out.history(:))), ...
    'The comparison did not complete an accepted full twentieth year.');
assert(isequal(size(leak.jointFreeze.history), ...
    [size(out.history,1),4*numel(out.parameters.leakage.jointAnglesRad)]), ...
    'Joint-freezing history is incomplete.');
assert(out.parameters.liningK==0 && leak.bodyVolume_m3_per_m==0 && ...
    abs(sum(leak.volumeByJoint_m3_per_m)-leak.totalVolume_m3_per_m)<1e-9, ...
    'The joint-only conservative partition failed.');
assert(numel(leak.dailyTotalVolume_m3_per_m)==365 && ...
    numel(leak.monthlyTotalVolume_m3_per_m)==12, ...
    'Daily or monthly integration is incomplete.');
end

function exportHistory(out,folder,name)
h=array2table(out.history,'VariableNames',out.historyColumns);
freeze=out.leakage.jointFreeze;
numberJoints=numel(out.parameters.leakage.jointAnglesRad);
for joint=1:numberJoints
    h.(sprintf('joint%d_factor',joint))=freeze.factor(:,joint);
    h.(sprintf('joint%d_minTemperature_C',joint))=freeze.minTemperature_C(:,joint);
    h.(sprintf('joint%d_maxTemperature_C',joint))=freeze.maxTemperature_C(:,joint);
    h.(sprintf('joint%d_frozenLengthFraction',joint))=freeze.frozenLengthFraction(:,joint);
    h.(sprintf('joint%d_conductance_m_per_s',joint))=freeze.conductance_m_per_s(:,joint);
end
writetable(h,fullfile(folder,[name '_history.csv']));
leak=out.leakage;
dayOfYear=leak.dailyBinEndDay(:)-19*365;
total_L_per_m=1000*leak.dailyTotalVolume_m3_per_m(:);
daily=table(dayOfYear,total_L_per_m);
for joint=1:numberJoints
    daily.(sprintf('joint%d_L_per_m',joint))=1000*leak.dailyVolumeByJoint_m3_per_m(:,joint);
end
writetable(daily,fullfile(folder,[name '_daily.csv']));
end
