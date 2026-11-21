function report=benchmark_TH_v6(checkpointFile,durationDays)
% BENCHMARK_TH_V6 Compare v5 full Newton and v6 adaptive modified Newton.
% Example:
%   report=benchmark_TH_v6('spinup60_output/active_year_FAILED_checkpoint.mat',0.001)
% Use a COPY of a checkpoint.  This function never overwrites it.
root=fileparts(mfilename('fullpath')); addpath(root);
if nargin<1 || isempty(checkpointFile)
    checkpointFile=fullfile(root,'spinup60_output','active_year_FAILED_checkpoint.mat');
end
if nargin<2, durationDays=.001; end
assert(isfile(checkpointFile),'Checkpoint file not found. Pass its full path.');
assert(durationDays>0 && durationDays<=1,'Benchmark duration must be in (0,1] day.');

fixedStep=durationDays*86400/8;
o=struct('meshFile',fullfile(root,'example_mesh_verified.mat'), ...
    'durationDays',durationDays,'dtInitial',fixedStep, ...
    'dtMax',fixedStep,'dtMin',fixedStep, ...
    'plot',false,'saveOutput',false,'verbose',false);
o.initial=struct('mode','checkpoint','checkpointFile',checkpointFile);
o.mechanics=struct('enable',false);
o.solver=struct('fd',1e-8,'maxIter',60,'tol',2e-8);

tic; a=TH_2D_conservative_v5(o); wallV5=toc;
o.solver.jacobianReuse=4; o.solver.jacobianRefreshRatio=.80;
tic; b=TH_2D_conservative_v6(o); wallV6=toc;

report=struct();
report.durationDays=durationDays; report.wallSecondsV5=wallV5;
report.wallSecondsV6=wallV6; report.speedup=wallV5/wallV6;
report.maxAbsDeltaT=max(abs(a.T-b.T));
report.maxAbsDeltaH=max(abs(a.h-b.h));
report.relativeLeakDifference=abs(sum(a.history(:,10).*a.history(:,2))- ...
    sum(b.history(:,10).*b.history(:,2)))/ ...
    max(abs(sum(a.history(:,10).*a.history(:,2))),realmin);
report.v5AcceptedSteps=size(a.history,1);
report.v5NonlinearIterations=sum(a.history(:,3));
report.v6Work=b.solverWork;
report.physicsRegressionPassed=report.maxAbsDeltaT<1e-5 && ...
    report.maxAbsDeltaH<1e-3 && report.relativeLeakDifference<1e-4;
assert(report.physicsRegressionPassed, ...
    'v5/v6 regression mismatch: inspect tolerances, fields and leakage before a long run.');
fprintf(['v5 %.3g s; v6 %.3g s; speed-up %.3gx; ' ...
    'max |dT| %.3g C; max |dh| %.3g m; leak rel. diff %.3g.\n'], ...
    wallV5,wallV6,report.speedup,report.maxAbsDeltaT,report.maxAbsDeltaH, ...
    report.relativeLeakDifference);
end
