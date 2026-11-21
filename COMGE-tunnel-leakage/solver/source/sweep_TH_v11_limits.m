function rows=sweep_TH_v11_limits(days)
%SWEEP_TH_V11_LIMITS Inspect trust-radius limits on year-20 continuation.
if nargin<1, days=30; end
root=fileparts(mfilename('fullpath'));
s=load(fullfile(root,'segment_v10_iteration_stress20_20260923','year_19.mat'),'out');
previous=s.out;
p=previous.parameters;
p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h,'startDay',19*365);
p.durationDays=days; p.outputEveryDays=days;
p.saveOutput=false; p.plot=false; p.verbose=false; p.mechanics.enable=false;
limits=[2 50;5 50;2 200;5 200;10 500];
rows=zeros(size(limits,1),9);
for k=1:size(limits,1)
    p.solver.stepT=limits(k,1); p.solver.stepH=limits(k,2);
    tic; out=TH_2D_segment_v11(p); elapsed=toc;
    rows(k,:)=[limits(k,:),elapsed,out.solverWork.acceptedSteps, ...
        out.solverWork.rejectedSteps,out.solverWork.jacobianBuilds, ...
        out.solverWork.nonlinearIterations, ...
        1000*out.leakage.totalVolume_m3_per_m,max(out.history(:,11))];
    fprintf('%g d: T=%.1f H=%.1f elapsed=%.2fs steps=%d rejects=%d J=%d iterations=%d leak=%.9f L/m\n', ...
        days,rows(k,1:8));
end
end
