function result=stress_test_TH_v9_target360_year20(outputFolder)
%STRESS_TEST_TH_V9_TARGET360_YEAR20 Full 20-year numerical A/B audit.
% One completed year is saved per call to the annual controller, so this
% function can be rerun after an interruption with the same source files.
% This is a numerical stress test, not a calibrated 20-year prediction.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'segment_v9_annular_target360_stress20_20260923');
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
checkpointFile=fullfile(outputFolder,'spinup_progress.mat');
marker=fullfile(outputFolder,'PAUSE_AFTER_YEAR');
c=struct('maxYears',20,'consecutiveYears',21,'outputFolder',outputFolder, ...
    'model',struct('solver',struct('dtTargetResidualEvals',360)));
result=[];
while true
    completed=0;
    if isfile(checkpointFile)
        saved=load(checkpointFile,'progress');
        completed=saved.progress.completedYears;
    end
    if completed>=20, break; end
    fid=fopen(marker,'w');
    assert(fid>0,'Cannot create annual pause marker.');
    fclose(fid);
    c.resume=completed>0;
    started=tic;
    result=run_TH_segment_spinup_v9(c);
    elapsed_s=toc(started);
    if isfile(marker), delete(marker); end
    if result.completedYears<=completed
        fprintf('Stress test stopped in year %d: %s\n',completed+1,result.status);
        break
    end
    assert(result.completedYears==completed+1,'Expected exactly one completed year.');
    year=result.completedYears;
    out=result.lastYear; %#ok<NASGU>
    assert(strcmp(out.version,'TH_v9_20260923_annular'));
    yearFile=fullfile(outputFolder,sprintf('year_%02d.mat',year));
    assert(~isfile(yearFile),'Refusing to overwrite a saved annual result.');
    save(yearFile,'out','-v7');
    h=out.history;
    phase=(h(:,1)-(year-1)*365)/(365/12);
    outputRemainder=abs(phase-round(phase))<1e-8;
    dt=h(:,2);
    row=table(year,1000*out.leakage.totalVolume_m3_per_m, ...
        out.solverWork.acceptedSteps,out.solverWork.rejectedSteps, ...
        min(dt),median(dt),sum(dt<.1),sum(dt<.01),sum(dt<.001), ...
        sum(dt<.1 & ~outputRemainder),sum(dt<.01 & ~outputRemainder), ...
        longestRun(dt<.1 & ~outputRemainder),max(h(:,3)),max(h(:,11)), ...
        elapsed_s,'VariableNames',{'Year','Leak_L_per_m','Accepted','Rejected', ...
        'MinDt_day','MedianDt_day','StepsBelow_0p1d','StepsBelow_0p01d', ...
        'StepsBelow_0p001d','InteriorBelow_0p1d','InteriorBelow_0p01d', ...
        'LongestInteriorRunBelow_0p1d','MaxNewtonIterations', ...
        'MaxScaledResidual','Elapsed_s'});
    statsFile=fullfile(outputFolder,'yearly_step_audit.mat');
    if isfile(statsFile)
        previous=load(statsFile,'audit'); audit=previous.audit;
        assert(~any(audit.Year==year),'Annual step audit already contains this year.');
        audit=[audit;row]; %#ok<AGROW>
    else
        audit=row;
    end
    save(statsFile,'audit','-v7');
    writetable(audit,fullfile(outputFolder,'yearly_step_audit.csv'));
    fprintf(['AUDIT year %d: %.6f L/m, %d accepted, %d rejected, ' ...
        'median dt %.5g day, min dt %.5g day, %d non-output steps <0.1 day\n'], ...
        year,row.Leak_L_per_m,row.Accepted,row.Rejected,row.MedianDt_day, ...
        row.MinDt_day,row.InteriorBelow_0p1d);
end
if isempty(result) && isfile(fullfile(outputFolder,'TH_segment_spinup_results.mat'))
    saved=load(fullfile(outputFolder,'TH_segment_spinup_results.mat'),'result');
    result=saved.result;
end
end

function n=longestRun(mask)
edges=diff([false;mask(:);false]);
starts=find(edges==1); ends=find(edges==-1);
if isempty(starts), n=0; else, n=max(ends-starts); end
end
