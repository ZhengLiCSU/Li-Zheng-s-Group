function result=stress_test_TH_v10_year20(outputFolder)
%STRESS_TEST_TH_V10_YEAR20 Reproducible 20-year numerical step audit.
% Rerun with the same output folder after interruption. Annual results are
% archived from the completed-year checkpoint before another year starts.
% This is a numerical stress test, not a calibrated 20-year prediction.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'segment_v10_iteration_stress20_20260923');
elseif ischar(outputFolder) && strcmp(outputFolder,'selftest')
    result=selftestRecovery();
    return
end
if ~isfolder(outputFolder), mkdir(outputFolder); end
checkpointFile=fullfile(outputFolder,'spinup_progress.mat');
marker=fullfile(outputFolder,'PAUSE_AFTER_YEAR');
ensureProvenance(root,outputFolder,checkpointFile);
c=struct('maxYears',20,'consecutiveYears',21,'outputFolder',outputFolder);
result=[];
while true
    completed=completedYears(checkpointFile);
    archiveCompleted(outputFolder,checkpointFile,completed,NaN);
    if completed>=20, break; end
    fid=fopen(marker,'w');
    assert(fid>0,'Cannot create annual pause marker.');
    fclose(fid);
    c.resume=completed>0;
    started=tic;
    try
        result=run_TH_segment_spinup_v10(c);
    catch ME
        removeMarker(marker);
        rethrow(ME);
    end
    elapsed_s=toc(started);
    removeMarker(marker);
    newCompleted=completedYears(checkpointFile);
    assert(newCompleted==result.completedYears, ...
        'Annual result and completed-year checkpoint disagree.');
    assert(newCompleted<=completed+1,'Expected at most one completed year.');
    if newCompleted>completed
        archiveCompleted(outputFolder,checkpointFile,newCompleted,elapsed_s);
    else
        fprintf('Stress test stopped in year %d: %s\n',completed+1,result.status);
        break
    end
end
if isempty(result) && isfile(fullfile(outputFolder,'TH_segment_spinup_results.mat'))
    saved=load(fullfile(outputFolder,'TH_segment_spinup_results.mat'),'result');
    result=saved.result;
end
end

function n=completedYears(checkpointFile)
n=0;
if isfile(checkpointFile)
    saved=load(checkpointFile,'progress');
    assert(isfield(saved,'progress') && isfield(saved.progress,'completedYears'), ...
        'Invalid completed-year checkpoint.');
    n=saved.progress.completedYears;
    assert(isscalar(n) && n==fix(n) && n>=0 && n<=20);
end
end

function archiveCompleted(folder,checkpointFile,completed,latestElapsed)
if completed==0, return; end
statsFile=fullfile(folder,'yearly_step_audit.mat');
csvFile=fullfile(folder,'yearly_step_audit.csv');
if isfile(statsFile)
    previous=load(statsFile,'audit');
    assert(isfield(previous,'audit') && istable(previous.audit), ...
        'Invalid annual audit file.');
    audit=previous.audit;
    assert(numel(unique(audit.Year))==height(audit) && ...
        all(audit.Year>=1 & audit.Year<=completed), ...
        'Annual audit has duplicate or out-of-range years.');
else
    audit=table();
end
for year=1:completed
    yearFile=fullfile(folder,sprintf('year_%02d.mat',year));
    recovered=false;
    elapsed_s=NaN;
    if isfile(yearFile)
        saved=load(yearFile,'out','elapsed_s');
        assert(isfield(saved,'out'),'Invalid archived year file.');
        out=saved.out;
        if isfield(saved,'elapsed_s'), elapsed_s=saved.elapsed_s; end
    else
        assert(year==completed, ...
            'An older year archive is missing; the latest checkpoint cannot recover it.');
        saved=load(checkpointFile,'progress');
        assert(saved.progress.completedYears==year, ...
            'Completed-year checkpoint changed during archival recovery.');
        out=saved.progress.lastYear;
        elapsed_s=latestElapsed;
        recovered=true;
        tempYear=[tempname(folder) '.mat'];
        save(tempYear,'out','elapsed_s','-v7');
        assert(~isfile(yearFile),'Refusing to overwrite a saved annual result.');
        movefile(tempYear,yearFile);
    end
    assert(strcmp(out.version,'TH_v10_20260923_iteration_dt') && ...
        abs(out.timeDay-year*365)<1e-7, ...
        'Annual archive has a wrong solver version or endpoint.');
    if isempty(audit) || ~any(audit.Year==year)
        row=auditRow(year,out,elapsed_s,recovered);
        if isempty(audit), audit=row; else, audit=[audit;row]; end %#ok<AGROW>
        fprintf(['AUDIT year %d: %.6f L/m, %d accepted, %d rejected, ' ...
            'median dt %.5g day, min dt %.5g day, %d non-output steps <0.1 day\n'], ...
            year,row.Leak_L_per_m,row.Accepted,row.Rejected,row.MedianDt_day, ...
            row.MinDt_day,row.InteriorBelow_0p1d);
    else
        k=find(audit.Year==year);
        assert(abs(audit.Leak_L_per_m(k)-1000*out.leakage.totalVolume_m3_per_m)<1e-8 ...
            && audit.Accepted(k)==out.solverWork.acceptedSteps, ...
            'Annual audit and year archive disagree.');
    end
end
audit=sortrows(audit,'Year');
tempStats=[tempname(folder) '.mat'];
save(tempStats,'audit','-v7');
movefile(tempStats,statsFile,'f');
tempCsv=[tempname(folder) '.csv'];
writetable(audit,tempCsv);
movefile(tempCsv,csvFile,'f');
end

function row=auditRow(year,out,elapsed_s,recovered)
h=out.history;
assert(size(h,2)==14 && ~isempty(h) && ...
    abs(sum(h(:,2))-365)<1e-6 && all(diff(h(:,1))>0), ...
    'Incomplete annual solver history.');
phase=(h(:,1)-(year-1)*365)/(365/12);
outputRemainder=abs(phase-round(phase))<1e-8;
dt=h(:,2);
row=table(year,1000*out.leakage.totalVolume_m3_per_m, ...
    out.solverWork.acceptedSteps,out.solverWork.rejectedSteps, ...
    min(dt),median(dt),sum(dt<.1),sum(dt<.01),sum(dt<.001), ...
    sum(dt<.1 & ~outputRemainder),sum(dt<.01 & ~outputRemainder), ...
    longestRun(dt<.1 & ~outputRemainder),max(h(:,3)),max(h(:,11)), ...
    elapsed_s,recovered, ...
    'VariableNames',{'Year','Leak_L_per_m','Accepted','Rejected', ...
    'MinDt_day','MedianDt_day','StepsBelow_0p1d','StepsBelow_0p01d', ...
    'StepsBelow_0p001d','InteriorBelow_0p1d','InteriorBelow_0p01d', ...
    'LongestInteriorRunBelow_0p1d','MaxNewtonIterations', ...
    'MaxScaledResidual','Elapsed_s','RecoveredFromCheckpoint'});
end

function ensureProvenance(root,folder,checkpointFile)
names={'TH_2D_segment_v10.m','run_TH_segment_spinup_v10.m', ...
    'stress_test_TH_v10_year20.m','summarize_segment_leakage.m', ...
    'example_mesh_repaired_interface.mat'};
sourceHashes=struct('name',names,'sha256',cell(size(names)));
for k=1:numel(names)
    sourceHashes(k).sha256=fileHash(fullfile(root,names{k}));
end
hashFile=fullfile(folder,'source_hashes.mat');
if isfile(hashFile)
    saved=load(hashFile,'sourceHashes');
    assert(isfield(saved,'sourceHashes') && isequal(saved.sourceHashes,sourceHashes), ...
        'Source or mesh changed during the 20-year stress test. Use a new output folder.');
else
    assert(~isfile(checkpointFile), ...
        'Cannot resume a long run without a source-hash record.');
    tempHash=[tempname(folder) '.mat'];
    save(tempHash,'sourceHashes','-v7');
    movefile(tempHash,hashFile);
end
jsonFile=fullfile(folder,'source_hashes.json');
fid=fopen(jsonFile,'w');
assert(fid>0,'Cannot write source-hash JSON file.');
clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(sourceHashes));
end

function digest=fileHash(file)
fid=fopen(file,'rb');
assert(fid>0,['Cannot read source file: ' file]);
clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes=fread(fid,inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes(:),'int8'));
raw=typecast(md.digest(),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
end

function removeMarker(marker)
if isfile(marker), delete(marker); end
end

function n=longestRun(mask)
edges=diff([false;mask(:);false]);
starts=find(edges==1); ends=find(edges==-1);
if isempty(starts), n=0; else, n=max(ends-starts); end
end

function checks=selftestRecovery()
folder=tempname; mkdir(folder);
clean=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
checkpointFile=fullfile(folder,'spinup_progress.mat');
ensureProvenance(fileparts(mfilename('fullpath')),folder,checkpointFile);
h=zeros(365,14);
h(:,1)=(1:365)'; h(:,2)=1; h(:,3)=1; h(:,11)=1e-9;
out=struct('version','TH_v10_20260923_iteration_dt','timeDay',365, ...
    'history',h,'leakage',struct('totalVolume_m3_per_m',.2), ...
    'solverWork',struct('acceptedSteps',365,'rejectedSteps',7));
progress=struct('completedYears',1,'lastYear',out); %#ok<NASGU>
save(checkpointFile,'progress','-v7');
ensureProvenance(fileparts(mfilename('fullpath')),folder,checkpointFile);
archiveCompleted(folder,checkpointFile,1,NaN);
yearFile=fullfile(folder,'year_01.mat');
assert(isfile(yearFile));
before=fileHash(yearFile);
auditFile=fullfile(folder,'yearly_step_audit.mat');
saved=load(auditFile,'audit');
assert(saved.audit.RecoveredFromCheckpoint && saved.audit.Rejected==7);
delete(auditFile);
delete(fullfile(folder,'yearly_step_audit.csv'));
archiveCompleted(folder,checkpointFile,1,NaN);
after=fileHash(yearFile);
assert(strcmp(before,after),'Recovery overwrote an existing annual result.');
saved=load(auditFile,'audit');
assert(~saved.audit.RecoveredFromCheckpoint && saved.audit.Rejected==7);
checks=struct('recoveredMissingYear',true,'recreatedMissingAudit',true, ...
    'preservedExistingYear',true,'checkedSourceHashes',true);
disp(checks);
end
