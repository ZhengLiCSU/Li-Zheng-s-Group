function rows=sweep_TH_v11_reuse(days)
%SWEEP_TH_V11_REUSE Inspect Jacobian refresh policy on fixed year-20 data.
if nargin<1, days=30; end
root=fileparts(mfilename('fullpath'));
s=load(fullfile(root,'segment_v10_iteration_stress20_20260923','year_19.mat'),'out');
previous=s.out;
p=previous.parameters;
p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h,'startDay',19*365);
p.durationDays=days; p.outputEveryDays=days;
p.saveOutput=false; p.plot=false; p.verbose=false; p.mechanics.enable=false;
policies=[4 .65;8 .8;12 .9;20 .95;100 .99];
rows=zeros(size(policies,1),9);
for k=1:size(policies,1)
    p.solver.jacobianReuse=policies(k,1);
    p.solver.jacobianRefreshRatio=policies(k,2);
    tic;
    out=TH_2D_segment_v11(p);
    elapsed=toc;
    rows(k,:)=[policies(k,:),elapsed,out.solverWork.acceptedSteps, ...
        out.solverWork.rejectedSteps,out.solverWork.jacobianBuilds, ...
        out.solverWork.nonlinearIterations, ...
        1000*out.leakage.totalVolume_m3_per_m,max(out.history(:,11))];
    fprintf('%g d: reuse=%g ratio=%.2f elapsed=%.2fs steps=%d rejects=%d J=%d iterations=%d leak=%.9f L/m\n', ...
        days,rows(k,1:8));
end
end
