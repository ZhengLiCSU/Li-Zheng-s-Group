function audit=benchmark_TH_v10_year20_isolation(outputFolder)
%BENCHMARK_TH_V10_YEAR20_ISOLATION Compare year 20 from the exact same state.
% v9 target-360 year 19 supplies T and h. The saved v9 year-20 trajectory
% supplies the reference and physical/solver parameters. Only the v10
% accepted-step controller differs; its legacy residual-evaluation target
% input is accepted for compatibility but does not control accepted steps.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v10_year20_isolation_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
sourceFolder=fullfile(root,'segment_v9_annular_target360_stress20_20260923');
s=load(fullfile(sourceFolder,'year_19.mat'),'out'); previous=s.out;
s=load(fullfile(sourceFolder,'year_20.mat'),'out'); reference=s.out;
assert(strcmp(previous.version,'TH_v9_20260923_annular'));
assert(strcmp(reference.version,'TH_v9_20260923_annular'));
assert(previous.timeDay==19*365 && reference.timeDay==20*365);
assert(isequal(previous.mesh,reference.mesh));
assert(isequal(previous.pAll,reference.pAll) && ...
    isequal(previous.soilNodes,reference.soilNodes));
assert(max(abs(reference.initialX-[previous.T;previous.h]))==0, ...
    'Saved reference year 20 did not start from saved year 19.');
resultFile=fullfile(outputFolder,'year_20_v10_same_initial.mat');
if isfile(resultFile)
    z=load(resultFile,'out'); out=z.out;
else
    p=reference.parameters;
    p.mesh=reference.mesh;
    p.initial=struct('mode','fields','T',previous.T,'h',previous.h, ...
        'startDay',19*365);
    p.durationDays=365;
    p.outputEveryDays=365/12;
    p.outputPrefix=fullfile(outputFolder,'active_year_20');
    p.saveOutput=true;
    p.plot=false;
    p.verbose=false;
    fprintf('v10 year-20 isolation: start day %.0f, v9 reference accepted=%d.\n', ...
        19*365,reference.solverWork.acceptedSteps);
    started=tic;
    out=TH_2D_segment_v10(p); %#ok<NASGU>
    elapsed_s=toc(started); %#ok<NASGU>
    save(resultFile,'out','elapsed_s','-v7');
end
assert(strcmp(out.version,'TH_v10_20260923_iteration_dt'));
assert(out.timeDay==20*365);
assert(max(abs(out.initialX-reference.initialX))==0);
assert(isequal(out.pAll,reference.pAll) && ...
    isequal(out.soilNodes,reference.soilNodes));
assert(abs(sum(out.history(:,2))-365)<1e-6);
assert(abs(sum(reference.history(:,2))-365)<1e-6);
old=stepRow(reference,'v9_target360');
new=stepRow(out,'v10_iteration');
stepAudit=[old;new];
annual=out.leakage.totalVolume_m3_per_m;
annualRef=reference.leakage.totalVolume_m3_per_m;
monthly=out.leakage.monthlyTotalVolume_m3_per_m(:);
monthlyRef=reference.leakage.monthlyTotalVolume_m3_per_m(:);
daily=out.leakage.dailyTotalVolume_m3_per_m(:);
dailyRef=reference.leakage.dailyTotalVolume_m3_per_m(:);
assert(numel(monthly)==12 && numel(monthlyRef)==12);
assert(numel(daily)==365 && numel(dailyRef)==365);
center=[out.mesh.xc out.mesh.yc]; xy=out.pAll;
roi=hypot(xy(:,1)-center(1),xy(:,2)-center(2))<=16;
roiH=roi(out.soilNodes);
dT=abs(out.T-reference.T); dH=abs(out.h-reference.h);
dW=abs(out.state.total-reference.state.total);
dI=abs(out.state.thetaI-reference.state.thetaI);
comparison=table(1000*annual,1000*annualRef, ...
    100*abs(annual-annualRef)/max(annualRef,eps), ...
    1000*max(abs(monthly-monthlyRef)), ...
    100*sum(abs(monthly-monthlyRef))/max(annualRef,eps), ...
    100*sum(abs(daily-dailyRef))/max(annualRef,eps), ...
    max(dT),max(dH),max(dT(roi)),max(dH(roiH)), ...
    max(dW(roiH)),max(dI(roiH)), ...
    'VariableNames',{'V10Leak_L_per_m','ReferenceLeak_L_per_m', ...
    'AnnualLeakRelativePercent','MaxMonthlyLeakDiff_L_per_m', ...
    'MonthlyLeakL1RelativePercent','DailyLeakL1RelativePercent', ...
    'MaxEndDeltaT_C','MaxEndDeltaH_m','RoiMaxEndDeltaT_C', ...
    'RoiMaxEndDeltaH_m','RoiMaxEndDeltaWater','RoiMaxEndDeltaIce'});
monthlyAudit=table((1:12)',1000*monthly,1000*monthlyRef, ...
    1000*(monthly-monthlyRef), ...
    100*abs(monthly-monthlyRef)./max(monthlyRef,1e-5), ...
    'VariableNames',{'Month','V10Leak_L_per_m','ReferenceLeak_L_per_m', ...
    'V10MinusReference_L_per_m','AbsRelativePercentWith0p01LFloor'});
writetable(stepAudit,fullfile(outputFolder,'year20_step_audit.csv'));
writetable(comparison,fullfile(outputFolder,'year20_comparison.csv'));
writetable(monthlyAudit,fullfile(outputFolder,'year20_monthly_audit.csv'));
audit=struct('stepAudit',stepAudit,'comparison',comparison, ...
    'monthlyAudit',monthlyAudit,'referenceVersion',reference.version, ...
    'candidateVersion',out.version);
save(fullfile(outputFolder,'year20_isolation_audit.mat'),'audit','-v7');
disp(stepAudit); disp(comparison);
end

function row=stepRow(out,name)
h=out.history; dt=h(:,2);
phase=(h(:,1)-19*365)/(365/12);
boundary=abs(phase-round(phase))<1e-8;
interior=dt(~boundary);
row=table(string(name),size(h,1),out.solverWork.rejectedSteps, ...
    min(dt),min(interior),median(dt), ...
    sum(dt<.1 & ~boundary),sum(dt<.01 & ~boundary), ...
    sum(dt<.001 & ~boundary),longestRun(dt<.1 & ~boundary), ...
    max(h(:,3)),max(h(:,11)), ...
    max(abs(h(:,8))),sum(abs(h(:,8))), ...
    max(abs(h(:,9))),sum(abs(h(:,9))), ...
    'VariableNames',{'Controller','Accepted','Rejected','MinDt_day', ...
    'MinInteriorDt_day','MedianDt_day','InteriorBelow_0p1d', ...
    'InteriorBelow_0p01d','InteriorBelow_0p001d', ...
    'LongestInteriorRunBelow_0p1d','MaxNewtonIterations', ...
    'MaxScaledResidual','MaxAbsWaterDefect_m2_per_m', ...
    'SumAbsWaterDefect_m2_per_m','MaxAbsEnergyDefect_J_per_m', ...
    'SumAbsEnergyDefect_J_per_m'});
end

function n=longestRun(mask)
edges=diff([false;mask(:);false]);
starts=find(edges==1); ends=find(edges==-1);
if isempty(starts), n=0; else, n=max(ends-starts); end
end
