function audit=benchmark_TH_v11_refresh_30d(outputFolder)
%BENCHMARK_TH_V11_REFRESH_30D Fixed fast4 year-20 start; vary only refresh.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v11_refresh_30d_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
maxNumCompThreads(16);
s=load(fullfile(root,'segment_v11_full20_fast4_20260923','year_19.mat'),'out');
previous=s.out;
p=previous.parameters;
p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h,'startDay',19*365);
p.durationDays=30; p.outputEveryDays=30;
p.saveOutput=false; p.plot=false; p.verbose=false; p.mechanics.enable=false;
assert(abs(p.dtMax/86400-4)<1e-12);
assert(p.solver.jacobianReuse==4 && abs(p.solver.jacobianRefreshRatio-.65)<1e-12);
ratios=[.65;.75;.85;.95;.65];
rows=zeros(numel(ratios),17);
outputs=cell(numel(ratios),1);
for k=1:numel(ratios)
    p.solver.jacobianRefreshRatio=ratios(k);
    started=tic;
    out=TH_2D_segment_v11(p);
    elapsed_s=toc(started);
    assert(isequal(out.initialX,[previous.T;previous.h]));
    assert(abs(out.timeDay-19*365-30)<1e-8 && abs(sum(out.history(:,2))-30)<1e-8);
    h=out.history;
    rows(k,:)=[ratios(k),elapsed_s, ...
        out.solverWork.acceptedSteps,out.solverWork.rejectedSteps, ...
        out.solverWork.nonlinearIterations,out.solverWork.jacobianBuilds, ...
        out.solverWork.linearSolves,out.solverWork.residualEvaluations, ...
        min(h(:,2)),median(h(:,2)),max(h(:,3)),max(h(:,11)), ...
        max(abs(h(:,8))),max(abs(h(:,9))), ...
        1000*out.leakage.totalVolume_m3_per_m, ...
        sum(abs(h(:,8))),sum(abs(h(:,9)))];
    outputs{k}=out;
    fprintf('refresh %.2f run %d: %.2fs, %d steps/%d rejects, Newton %d, J %d, solves %d, residuals %d, leak %.9f L/m\n', ...
        ratios(k),k,elapsed_s,rows(k,3:8),rows(k,15));
    save(fullfile(outputFolder,sprintf('ratio_%02d_run_%d.mat',round(100*ratios(k)),k)), ...
        'out','elapsed_s','-v7');
end
base=outputs{1};
deltaT=zeros(numel(ratios),1); deltaH=deltaT; deltaLeak=deltaT;
for k=1:numel(ratios)
    out=outputs{k};
    deltaT(k)=max(abs(out.T-base.T));
    deltaH(k)=max(abs(out.h-base.h));
    deltaLeak(k)=100*abs(out.leakage.totalVolume_m3_per_m- ...
        base.leakage.totalVolume_m3_per_m)/max(base.leakage.totalVolume_m3_per_m,eps);
end
audit=array2table([rows,deltaT,deltaH,deltaLeak], ...
    'VariableNames',{'RefreshRatio','Elapsed_s','Accepted','Rejected', ...
    'NewtonIterations','JacobianBuilds','LinearSolves','ResidualEvaluations', ...
    'MinDt_day','MedianDt_day','MaxStepNewton','MaxScaledResidual', ...
    'MaxAbsWaterDefect_m2_per_m','MaxAbsEnergyDefect_J_per_m', ...
    'Leak_L_per_m','SumAbsWaterDefect_m2_per_m','SumAbsEnergyDefect_J_per_m', ...
    'MaxEndDeltaT_vsDefault_C','MaxEndDeltaH_vsDefault_m', ...
    'LeakRelative_vsDefault_percent'});
writetable(audit,fullfile(outputFolder,'refresh_30d_audit.csv'));
save(fullfile(outputFolder,'refresh_30d_audit.mat'),'audit','-v7');
disp(audit);
end
