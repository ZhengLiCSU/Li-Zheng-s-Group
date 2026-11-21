function report=benchmark_TH_v10_timestep_ref(outputFolder)
%BENCHMARK_TH_V10_TIMESTEP_REF Check year-4 temporal resolution at 0.5 day.
% Start from the same saved v9 year-3 state as the v10 A/B benchmark.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v10_ab_from_v9_year3_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
s=load(fullfile(root,'segment_v9_annular_stress20_20260923','year_03.mat'),'out');
previous=s.out;
assert(strcmp(previous.version,'TH_v9_20260923_annular'));
v=load(fullfile(outputFolder,'year_04_v10.mat'),'out');
base=v.out;
assert(strcmp(base.version,'TH_v10_20260923_iteration_dt'));
assert(max(abs(base.initialX-[previous.T;previous.h]))<1e-12);
resultFile=fullfile(outputFolder,'year_04_v10_dtmax_halfday.mat');
if isfile(resultFile)
    z=load(resultFile,'half'); half=z.half;
else
    p=previous.parameters;
    p.mesh=previous.mesh;
    p.initial=struct('mode','fields','T',previous.T,'h',previous.h, ...
        'startDay',1095);
    p.durationDays=365;
    p.outputEveryDays=365/12;
    p.dtInitial=.1*86400;
    p.dtMax=.5*86400;
    p.saveOutput=true;
    p.plot=false;
    p.verbose=false;
    p.mechanics.enable=false;
    p.outputPrefix=fullfile(outputFolder,'active_year_04_dtmax_halfday');
    fprintf('v10 temporal reference: starting year 4 with dtMax=0.5 day.\n');
    started=tic;
    half=TH_2D_segment_v10(p); %#ok<NASGU>
    elapsed_s=toc(started); %#ok<NASGU>
    save(resultFile,'half','elapsed_s','-v7');
end
assert(strcmp(half.version,'TH_v10_20260923_iteration_dt'));
assert(abs(half.timeDay-1460)<1e-7);
assert(isequal(base.pAll,half.pAll) && isequal(base.soilNodes,half.soilNodes));
xy=base.pAll; center=[base.mesh.xc base.mesh.yc];
roi=hypot(xy(:,1)-center(1),xy(:,2)-center(2))<=16;
roiH=roi(base.soilNodes);
dT=abs(base.T-half.T); dH=abs(base.h-half.h);
monthly=base.leakage.monthlyTotalVolume_m3_per_m;
monthlyHalf=half.leakage.monthlyTotalVolume_m3_per_m;
daily=base.leakage.dailyTotalVolume_m3_per_m;
dailyHalf=half.leakage.dailyTotalVolume_m3_per_m;
annual=base.leakage.totalVolume_m3_per_m;
annualHalf=half.leakage.totalVolume_m3_per_m;
assert(numel(monthly)==12 && numel(monthlyHalf)==12);
assert(numel(daily)==365 && numel(dailyHalf)==365);
monthlyBase_L=1000*monthly(:);
monthlyHalf_L=1000*monthlyHalf(:);
monthlyAudit=table((1:12)',monthlyBase_L,monthlyHalf_L, ...
    monthlyBase_L-monthlyHalf_L, ...
    100*abs(monthlyBase_L-monthlyHalf_L)./max(monthlyHalf_L,.01), ...
    'VariableNames',{'Month','Base_L_per_m','HalfDay_L_per_m', ...
    'BaseMinusHalfDay_L_per_m','AbsRelativePercentWith0p01LFloor'});
writetable(monthlyAudit,fullfile(outputFolder,'v10_halfday_monthly_audit.csv'));
dt=half.history(:,2);
phase=(half.history(:,1)-1095)/(365/12);
outputRemainder=abs(phase-round(phase))<1e-8;
report=table(base.solverWork.acceptedSteps,half.solverWork.acceptedSteps, ...
    base.solverWork.rejectedSteps,half.solverWork.rejectedSteps, ...
    1000*annual,1000*annualHalf, ...
    100*abs(annual-annualHalf)/max(annualHalf,eps), ...
    1000*max(abs(monthly(:)-monthlyHalf(:))), ...
    100*sum(abs(daily(:)-dailyHalf(:)))/max(annualHalf,eps), ...
    max(dT),max(dH),max(dT(roi)),max(dH(roiH)), ...
    min(dt),sum(dt<.1 & ~outputRemainder), ...
    max(abs(half.history(:,8))),max(abs(half.history(:,9))), ...
    'VariableNames',{'BaseAccepted','HalfDayAccepted','BaseRejected', ...
    'HalfDayRejected','BaseLeak_L_per_m','HalfDayLeak_L_per_m', ...
    'AnnualLeakRelativePercent','MaxMonthlyLeakDiff_L_per_m', ...
    'DailyLeakL1RelativePercent','MaxEndDeltaT_C','MaxEndDeltaH_m', ...
    'RoiMaxEndDeltaT_C','RoiMaxEndDeltaH_m','HalfDayMinDt_day', ...
    'HalfDayInteriorBelow_0p1d','HalfDayMaxAbsWaterDefect_m2_per_m', ...
    'HalfDayMaxAbsEnergyDefect_J_per_m'});
writetable(report,fullfile(outputFolder,'v10_halfday_timestep_audit.csv'));
save(fullfile(outputFolder,'v10_halfday_timestep_audit.mat'),'report','-v7');
disp(report);
end
