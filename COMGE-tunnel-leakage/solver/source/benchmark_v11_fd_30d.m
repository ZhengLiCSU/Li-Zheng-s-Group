function audit=benchmark_v11_fd_30d(outputFolder)
%BENCHMARK_V11_FD_30D Same-start tangent-perturbation sweep; residual unchanged.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v11_fd_sweep_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
maxNumCompThreads(16);
s=load(fullfile(root,'segment_v11_full20_fast4_20260923','year_19.mat'),'out');
previous=s.out; p=previous.parameters; p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h,'startDay',19*365);
p.durationDays=30; p.outputEveryDays=30;
p.saveOutput=false; p.plot=false; p.verbose=false; p.mechanics.enable=false;
steps=[1e-8;1e-7;1e-6;1e-5;1e-8];
rows=zeros(numel(steps),14); results=cell(numel(steps),1);
for k=1:numel(steps)
    p.solver.fd=steps(k);
    started=tic; out=TH_2D_segment_v11(p); elapsed=toc(started);
    results{k}=out; h=out.history;
    rows(k,:)=[k,steps(k),elapsed,out.solverWork.acceptedSteps, ...
        out.solverWork.rejectedSteps,out.solverWork.nonlinearIterations, ...
        out.solverWork.jacobianBuilds,out.solverWork.linearSolves, ...
        out.solverWork.residualEvaluations, ...
        1000*out.leakage.totalVolume_m3_per_m,max(h(:,11)), ...
        max(abs(h(:,8))),max(abs(h(:,9))),min(h(:,2))];
    fprintf('fd %.0e run%d: %.3fs, Newton %d J %d residual %d Q %.9f L/m\n', ...
        steps(k),k,elapsed,out.solverWork.nonlinearIterations, ...
        out.solverWork.jacobianBuilds,out.solverWork.residualEvaluations,rows(k,10));
end
audit=array2table(rows,'VariableNames',{'Run','Fd','Elapsed_s','Accepted', ...
    'Rejected','Newton','Jacobians','LinearSolves','Residuals', ...
    'Leak_L_per_m','MaxScaledResidual','MaxWaterDefect_m2_per_m', ...
    'MaxEnergyDefect_J_per_m','MinDt_day'});
base=results{1};
delta=zeros(numel(steps),3);
for k=1:numel(steps)
    out=results{k};
    delta(k,:)=[max(abs(out.T-base.T)),max(abs(out.h-base.h)), ...
        100*abs(out.leakage.totalVolume_m3_per_m- ...
        base.leakage.totalVolume_m3_per_m)/base.leakage.totalVolume_m3_per_m];
end
audit.EndDeltaT_C=delta(:,1); audit.EndDeltaH_m=delta(:,2);
audit.LeakRelativeError_percent=delta(:,3);
writetable(audit,fullfile(outputFolder,'fd_30d_audit.csv'));
save(fullfile(outputFolder,'fd_30d_results.mat'),'audit','-v7');
disp(audit);
end
