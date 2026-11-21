function report=benchmark_TH_v10_year20_dtmax(dtMaxDays)
%BENCHMARK_TH_V10_YEAR20_DTMAX Time coarser year-20 stepping from one state.
if nargin<1, dtMaxDays=2; end
assert(isscalar(dtMaxDays) && dtMaxDays>1 && isfinite(dtMaxDays));
root=fileparts(mfilename('fullpath'));
source=fullfile(root,'segment_v10_iteration_stress20_20260923');
s=load(fullfile(source,'year_19.mat'),'out'); previous=s.out;
s=load(fullfile(source,'year_20.mat'),'out'); base=s.out;
assert(isequal(base.initialX,[previous.T;previous.h]), ...
    'The reference year does not start from the saved year-19 state.');
folder=fullfile(root,'v10_year20_dtmax_sweep');
if ~isfolder(folder), mkdir(folder); end
name=strrep(sprintf('%.5g',dtMaxDays),'.','p');
file=fullfile(folder,['year20_dtmax_' name 'day.mat']);
if isfile(file)
    s=load(file,'test','elapsed_s'); test=s.test; elapsed_s=s.elapsed_s;
else
    p=previous.parameters;
    p.mesh=previous.mesh;
    p.initial=struct('mode','fields','T',previous.T, ...
        'h',previous.h,'startDay',19*365);
    p.durationDays=365; p.outputEveryDays=365/12;
    p.dtInitial=.1*86400; p.dtMax=dtMaxDays*86400;
    p.saveOutput=false; p.plot=false; p.verbose=false;
    p.mechanics.enable=false;
    started=tic;
    test=TH_2D_segment_v10(p); %#ok<NASGU>
    elapsed_s=toc(started); %#ok<NASGU>
    save(file,'test','elapsed_s','-v7');
end
assert(abs(test.timeDay-20*365)<1e-7 && ...
    isequal(test.initialX,base.initialX));
annualBase=1000*base.leakage.totalVolume_m3_per_m;
annualTest=1000*test.leakage.totalVolume_m3_per_m;
monthBase=1000*base.leakage.monthlyTotalVolume_m3_per_m(:);
monthTest=1000*test.leakage.monthlyTotalVolume_m3_per_m(:);
dailyBase=base.leakage.dailyTotalVolume_m3_per_m(:);
dailyTest=test.leakage.dailyTotalVolume_m3_per_m(:);
report=table(dtMaxDays,elapsed_s,size(test.history,1), ...
    test.solverWork.rejectedSteps,median(test.history(:,2)), ...
    annualBase,annualTest,100*abs(annualTest-annualBase)/annualBase, ...
    max(abs(monthTest-monthBase)), ...
    100*sum(abs(dailyTest-dailyBase))/base.leakage.totalVolume_m3_per_m, ...
    max(abs(test.T-base.T)),max(abs(test.h-base.h)), ...
    max(abs(test.history(:,8))),max(abs(test.history(:,9))), ...
    'VariableNames',{'DtMax_day','Elapsed_s','Accepted','Rejected', ...
    'MedianDt_day','BaseLeak_L_per_m','TestLeak_L_per_m', ...
    'AnnualLeakRelativePercent','MaxMonthlyLeakDiff_L_per_m', ...
    'DailyLeakL1RelativePercent','MaxEndDeltaT_C','MaxEndDeltaH_m', ...
    'MaxAbsWaterDefect_m2_per_m','MaxAbsEnergyDefect_J_per_m'});
writetable(report,fullfile(folder,['audit_' name 'day.csv']));
disp(report);
end
