function audit=benchmark_TH_v10_year4_5(outputFolder)
%BENCHMARK_TH_V10_YEAR4_5 Compare accepted-step control from a common year-3 state.
% Both v10 years start from the preceding v10 result; the reference is the
% independently saved v9 target-360 year-4/5 continuation. This benchmark
% changes time-step selection only and does not assert a physical calibration.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v10_ab_from_v9_year3_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
sourceFolder=fullfile(root,'segment_v9_annular_stress20_20260923');
source=load(fullfile(sourceFolder,'year_03.mat'),'out');
assert(strcmp(source.out.version,'TH_v9_20260923_annular'));
references=cell(2,1);
r=load(fullfile(sourceFolder,'year_04_target360.mat'),'alt');
references{1}=r.alt;
r=load(fullfile(sourceFolder,'year_05_target360.mat'),'alt5');
references{2}=r.alt5;
assert(all(cellfun(@(x) strcmp(x.version,'TH_v9_20260923_annular'),references)));
previous=source.out;
rows=cell(2,1);
for year=4:5
    ix=year-3;
    resultFile=fullfile(outputFolder,sprintf('year_%02d_v10.mat',year));
    if isfile(resultFile)
        z=load(resultFile,'out'); out=z.out;
        assert(strcmp(out.version,'TH_v10_20260923_iteration_dt'));
        assert(max(abs(out.initialX-[previous.T;previous.h]))<1e-12);
    else
        p=previous.parameters;
        p.mesh=previous.mesh;
        p.initial=struct('mode','fields','T',previous.T,'h',previous.h, ...
            'startDay',(year-1)*365);
        p.durationDays=365;
        p.outputEveryDays=365/12;
        p.dtInitial=.1*86400;
        p.dtMax=86400;
        p.saveOutput=true;
        p.plot=false;
        p.verbose=false;
        p.mechanics.enable=false;
        p.outputPrefix=fullfile(outputFolder,sprintf('active_year_%02d',year));
        fprintf('v10 A/B: starting year %d from day %.0f.\n',year,(year-1)*365);
        started=tic;
        out=TH_2D_segment_v10(p); %#ok<NASGU>
        elapsed_s=toc(started); %#ok<NASGU>
        assert(strcmp(out.version,'TH_v10_20260923_iteration_dt'));
        assert(abs(out.timeDay-year*365)<1e-7);
        assert(abs(sum(out.history(:,2))-365)<1e-6);
        save(resultFile,'out','elapsed_s','-v7');
    end
    ref=references{ix};
    assert(isequal(out.pAll,ref.pAll) && isequal(out.soilNodes,ref.soilNodes), ...
        'A/B mesh node ordering differs.');
    h=out.history; dt=h(:,2);
    phase=(h(:,1)-(year-1)*365)/(365/12);
    outputRemainder=abs(phase-round(phase))<1e-8;
    interior=dt(~outputRemainder);
    if isempty(interior), minInterior=NaN; else, minInterior=min(interior); end
    m=out.leakage.monthlyTotalVolume_m3_per_m;
    mr=ref.leakage.monthlyTotalVolume_m3_per_m;
    assert(numel(m)==12 && numel(mr)==12);
    annual=out.leakage.totalVolume_m3_per_m;
    annualRef=ref.leakage.totalVolume_m3_per_m;
    daily=out.leakage.dailyTotalVolume_m3_per_m;
    dailyRef=ref.leakage.dailyTotalVolume_m3_per_m;
    assert(numel(daily)==365 && numel(dailyRef)==365);
    xy=out.pAll; center=[out.mesh.xc out.mesh.yc];
    roi=hypot(xy(:,1)-center(1),xy(:,2)-center(2))<=16;
    roiH=roi(out.soilNodes);
    dT=abs(out.T-ref.T); dH=abs(out.h-ref.h);
    dW=abs(out.state.total-ref.state.total);
    dI=abs(out.state.thetaI-ref.state.thetaI);
    yearResult=table(year,1000*annual,1000*annualRef, ...
        100*abs(annual-annualRef)/max(annualRef,eps), ...
        1000*max(abs(m(:)-mr(:))), ...
        100*sum(abs(daily(:)-dailyRef(:)))/max(annualRef,eps), ...
        size(h,1),out.solverWork.rejectedSteps,min(dt),minInterior,median(dt), ...
        sum(dt<.1 & ~outputRemainder),sum(dt<.01 & ~outputRemainder), ...
        longestRun(dt<.1 & ~outputRemainder),max(h(:,3)), ...
        max(abs(h(:,8))),sum(abs(h(:,8))), ...
        max(abs(h(:,9))),sum(abs(h(:,9))), ...
        max(dT),max(dH),max(dT(roi)),max(dH(roiH)), ...
        max(dW(roiH)),max(dI(roiH)), ...
        'VariableNames',{'Year','Leak_L_per_m','ReferenceLeak_L_per_m', ...
        'AnnualLeakRelativePercent','MaxMonthlyLeakDiff_L_per_m', ...
        'DailyLeakL1RelativePercent','Accepted','Rejected','MinDt_day', ...
        'MinInteriorDt_day','MedianDt_day','InteriorBelow_0p1d', ...
        'InteriorBelow_0p01d','LongestInteriorRunBelow_0p1d', ...
        'MaxNewtonIterations','MaxAbsWaterDefect_m2_per_m', ...
        'SumAbsWaterDefect_m2_per_m','MaxAbsEnergyDefect_J_per_m', ...
        'SumAbsEnergyDefect_J_per_m','MaxEndDeltaT_C', ...
        'MaxEndDeltaH_m','RoiMaxEndDeltaT_C','RoiMaxEndDeltaH_m', ...
        'RoiMaxEndDeltaWater','RoiMaxEndDeltaIce'});
    rows{ix}=yearResult;
    audit=vertcat(rows{1:ix}); %#ok<NASGU>
    writetable(audit,fullfile(outputFolder,'v10_year4_5_ab_audit.csv'));
    save(fullfile(outputFolder,'v10_year4_5_ab_audit.mat'),'audit','-v7');
    fprintf('v10 A/B year %d: %d accepted, %d rejected, %d interior dt<0.1 d, ', ...
        year,yearResult.Accepted,yearResult.Rejected,yearResult.InteriorBelow_0p1d);
    fprintf('leak %.6f L/m, difference %.4f%%, ROI max dT %.4g C, dH %.4g m.\n', ...
        yearResult.Leak_L_per_m,yearResult.AnnualLeakRelativePercent, ...
        yearResult.RoiMaxEndDeltaT_C,yearResult.RoiMaxEndDeltaH_m);
    previous=out;
end
end

function n=longestRun(mask)
edges=diff([false;mask(:);false]);
starts=find(edges==1); ends=find(edges==-1);
if isempty(starts), n=0; else, n=max(ends-starts); end
end
