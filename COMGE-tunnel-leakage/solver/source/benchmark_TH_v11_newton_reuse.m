function audit=benchmark_TH_v11_newton_reuse(outputFolder)
%BENCHMARK_TH_V11_NEWTON_REUSE Same-state 30-day Jacobian refresh study.
% The v11 solver is not modified. Runs are serial, with identical model,
% initial fields, time span, mesh, tolerances, and timestep policy.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v11_newton_reuse_30d_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
s=load(fullfile(root,'segment_v10_iteration_stress20_20260923','year_19.mat'),'out');
previous=s.out;
p=previous.parameters;
p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h, ...
    'startDay',19*365);
p.durationDays=30; p.outputEveryDays=30;
p.saveOutput=false; p.plot=false; p.verbose=false;
p.mechanics.enable=false;
reuseValues=[0;1;2;4];
rows=zeros(numel(reuseValues),15);
results=cell(numel(reuseValues),1);
for k=1:numel(reuseValues)
    p.solver.jacobianReuse=reuseValues(k);
    started=tic;
    out=TH_2D_segment_v11(p);
    elapsed_s=toc(started);
    assert(max(abs(out.initialX-[previous.T;previous.h]))==0);
    assert(abs(out.timeDay-(19*365+30))<1e-8);
    assert(abs(sum(out.history(:,2))-30)<1e-8);
    h=out.history;
    rows(k,:)=[reuseValues(k),elapsed_s,size(h,1), ...
        out.solverWork.rejectedSteps,out.solverWork.nonlinearIterations, ...
        out.solverWork.jacobianBuilds,out.solverWork.linearSolves, ...
        out.solverWork.residualEvaluations, ...
        min(h(:,2)),median(h(:,2)), ...
        1000*out.leakage.totalVolume_m3_per_m,max(h(:,11)), ...
        max(abs(h(:,8))),max(abs(h(:,9))),sum(abs(h(:,8)))];
    results{k}=out;
    fprintf('reuse=%d: %.2fs, %d accepted/%d rejected, %d Newton, %d J/LU, %d residual, leak %.9f L/m\n', ...
        reuseValues(k),elapsed_s,rows(k,3:6),rows(k,8),rows(k,11));
    save(fullfile(outputFolder,sprintf('reuse_%d.mat',reuseValues(k))), ...
        'out','elapsed_s','-v7');
end
base=results{end};
leakBase=base.leakage.totalVolume_m3_per_m;
deltaT=zeros(numel(reuseValues),1); deltaH=deltaT; deltaLeak=deltaT;
for k=1:numel(reuseValues)
    out=results{k};
    deltaT(k)=max(abs(out.T-base.T));
    deltaH(k)=max(abs(out.h-base.h));
    deltaLeak(k)=100*abs(out.leakage.totalVolume_m3_per_m-leakBase)/max(leakBase,eps);
end
audit=array2table([rows,deltaT,deltaH,deltaLeak], ...
    'VariableNames',{'JacobianReuse','Elapsed_s','Accepted','Rejected', ...
    'NewtonIterations','JacobianBuilds','LinearSolves', ...
    'ResidualEvaluations','MinDt_day','MedianDt_day', ...
    'Leak_L_per_m','MaxScaledResidual','MaxWaterDefect_m2_per_m', ...
    'MaxEnergyDefect_J_per_m','SumAbsWaterDefect_m2_per_m', ...
    'MaxDeltaT_vsReuse4_C','MaxDeltaH_vsReuse4_m', ...
    'LeakRelative_vsReuse4_percent'});
writetable(audit,fullfile(outputFolder,'newton_reuse_audit.csv'));
save(fullfile(outputFolder,'newton_reuse_audit.mat'),'audit','-v7');
disp(audit);
end
