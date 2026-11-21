function audit=benchmark_v11_nonmonotone_abba(outputFolder,mode)
%BENCHMARK_V11_NONMONOTONE_ABBA Balanced repeated 30-day same-start timing.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'v11_nonmonotone_abba_20260923');
end
if nargin<2, mode='nonmonotone'; end
if ~isfolder(outputFolder), mkdir(outputFolder); end
maxNumCompThreads(16);
s=load(fullfile(root,'segment_v11_full20_fast4_20260923','year_19.mat'),'out');
previous=s.out;
p=previous.parameters; p.mesh=previous.mesh;
p.initial=struct('mode','fields','T',previous.T,'h',previous.h,'startDay',19*365);
p.durationDays=30; p.outputEveryDays=30;
p.saveOutput=false; p.plot=false; p.verbose=false; p.mechanics.enable=false;
pCandidate=p;
if strcmp(mode,'inf_merit')
    pCandidate.solver.nonmonotoneWindow=0;
    pCandidate.solver.meritNorm='inf';
else
    assert(strcmp(mode,'nonmonotone'));
end
methods={'base';'candidate';'candidate';'base'};
rows=zeros(4,13); results=cell(4,1);
for k=1:4
    started=tic;
    if strcmp(methods{k},'base')
        results{k}=TH_2D_segment_v11(p);
    else
        results{k}=TH_2D_segment_v11_nonlinear(pCandidate);
    end
    elapsed=toc(started); out=results{k}; h=out.history;
    rows(k,:)=[k,elapsed,out.solverWork.acceptedSteps, ...
        out.solverWork.rejectedSteps,out.solverWork.nonlinearIterations, ...
        out.solverWork.jacobianBuilds,out.solverWork.linearSolves, ...
        out.solverWork.residualEvaluations, ...
        1000*out.leakage.totalVolume_m3_per_m,max(h(:,11)), ...
        max(abs(h(:,8))),max(abs(h(:,9))),min(h(:,2))];
    fprintf('%d %s %.3fs; Newton %d J %d residual %d\n', ...
        k,methods{k},elapsed,out.solverWork.nonlinearIterations, ...
        out.solverWork.jacobianBuilds,out.solverWork.residualEvaluations);
end
audit=array2table(rows,'VariableNames',{'Run','Elapsed_s','Accepted', ...
    'Rejected','Newton','Jacobians','LinearSolves','Residuals', ...
    'Leak_L_per_m','MaxScaledResidual','MaxWaterDefect_m2_per_m', ...
    'MaxEnergyDefect_J_per_m','MinDt_day'});
audit.Method=methods; audit=movevars(audit,'Method','After','Run');
writetable(audit,fullfile(outputFolder,'abba_30d_audit.csv'));
b=results{1}; c=results{2};
if size(b.history,1)==size(c.history,1)
    dtMaxDiff=max(abs(b.history(:,2)-c.history(:,2)));
    timeMaxDiff=max(abs(b.history(:,1)-c.history(:,1)));
else
    dtMaxDiff=NaN; timeMaxDiff=NaN;
end
comparison=struct('baseMeanTime_s',mean(rows([1 4],2)), ...
    'candidateMeanTime_s',mean(rows([2 3],2)), ...
    'speedup',mean(rows([1 4],2))/mean(rows([2 3],2)), ...
    'maxStepDtDifference_day',dtMaxDiff,'maxStepEndTimeDifference_day',timeMaxDiff, ...
    'endDeltaT_C',max(abs(b.T-c.T)), ...
    'endDeltaH_m',max(abs(b.h-c.h)), ...
    'leakRelativeError_percent',100*abs(b.leakage.totalVolume_m3_per_m- ...
        c.leakage.totalVolume_m3_per_m)/b.leakage.totalVolume_m3_per_m);
save(fullfile(outputFolder,'abba_30d_results.mat'),'b','c','audit','comparison','-v7');
disp(audit); disp(comparison);
end
