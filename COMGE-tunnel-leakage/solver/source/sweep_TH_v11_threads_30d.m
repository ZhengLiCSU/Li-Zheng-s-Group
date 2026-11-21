function audit=sweep_TH_v11_threads_30d()
%SWEEP_TH_V11_THREADS_30D Compare thread settings on identical age-20 fields.
root=fileparts(mfilename('fullpath'));
reference=load(fullfile(root,'segment_v11_full20_fast4_20260923', ...
    'year_19.mat'),'out');
previous=reference.out;
p=previous.parameters;
p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h, ...
    'startDay',19*365);
p.durationDays=30;
p.outputEveryDays=30;
p.saveOutput=false;
p.plot=false;
p.verbose=false;
p.mechanics.enable=false;
order=[16 8 4 1 1 4 8 16];
data=zeros(numel(order),8);
for k=1:numel(order)
    maxNumCompThreads(order(k));
    started=tic;
    out=TH_2D_segment_v11(p);
    elapsed=toc(started);
    data(k,:)=[k order(k) elapsed out.solverWork.acceptedSteps ...
        out.solverWork.nonlinearIterations out.solverWork.jacobianBuilds ...
        out.solverWork.residualEvaluations ...
        1000*out.leakage.totalVolume_m3_per_m];
    fprintf('threads=%d run=%d elapsed=%.3f s Q=%.12g L/m\n', ...
        order(k),k,elapsed,data(k,8));
end
audit=array2table(data,'VariableNames',{'Run','Threads','Elapsed_s', ...
    'Accepted','Newton','Jacobians','Residuals','Leak_L_per_m'});
assert(max(abs(audit.Leak_L_per_m-audit.Leak_L_per_m(1)))<1e-8);
folder=fullfile(root,'v11_threads_30d_20260923');
if ~isfolder(folder), mkdir(folder); end
writetable(audit,fullfile(folder,'threads_audit.csv'));
end
