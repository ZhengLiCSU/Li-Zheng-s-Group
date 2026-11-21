function rows=profile_TH_v10_short(days)
%PROFILE_TH_V10_SHORT Profile a representative year-20 continuation.
if nargin<1, days=30; end
root=fileparts(mfilename('fullpath'));
s=load(fullfile(root,'segment_v10_iteration_stress20_20260923','year_19.mat'),'out');
previous=s.out;
p=previous.parameters;
p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h,'startDay',19*365);
p.durationDays=days;
p.outputEveryDays=days;
p.saveOutput=false; p.plot=false; p.verbose=false; p.mechanics.enable=false;
profile clear; profile on;
started=tic;
out=TH_2D_segment_v10(p);
elapsed=toc(started);
profile off;
q=profile('info');
f=q.FunctionTable;
[~,order]=sort([f.TotalTime],'descend');
n=min(numel(order),25);
rows=cell(n,4);
for i=1:n
    k=order(i);
    rows(i,:)={f(k).FunctionName,f(k).NumCalls,f(k).TotalTime,f(k).TotalTime/max(elapsed,eps)};
end
fprintf('days=%g elapsed=%.3f s, accepted=%d, Jacobians=%d, residual evaluations=%d\n', ...
    days,elapsed,out.solverWork.acceptedSteps,out.solverWork.jacobianBuilds,out.solverWork.residualEvaluations);
disp(cell2table(rows,'VariableNames',{'Function','Calls','Seconds','Fraction'}));
end
