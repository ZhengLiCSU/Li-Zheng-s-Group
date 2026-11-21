function result=run_TH_segment_spinup_v12(user)
% RUN_TH_SEGMENT_SPINUP_V12 Annual TH continuation with joint ice blockage.
% Run result=run_TH_segment_spinup_v12; resume with struct('resume',true).
% TH_2D_segment_v12.m and example_mesh_repaired_interface.mat must be on this path.
% All thresholds below are numerical research criteria, not design standards.
if nargin<1, user=struct(); end
root=fileparts(mfilename('fullpath')); addpath(root);
c=struct('maxYears',60,'consecutiveYears',3,'productionYears',1, ...
 'tolT',.02,'tolH',.05,'tolWater',1e-4,'tolIce',1e-4, ...
 'roiTolT',.05,'roiTolH',.10,'roiTolWater',5e-4,'roiTolIce',5e-4,'roiRadius',16, ...
 'tolLeakRelative',.005,'tolLeakAbsolute',1e-6, ...
 'resume',false,'outputFolder',fullfile(root,'segment_spinup_v12_output'), ...
 'model',struct());
names=fieldnames(user);
for k=1:numel(names)
 assert(isfield(c,names{k}),['Unknown option: ' names{k}]); c.(names{k})=user.(names{k});
end
assert(isscalar(c.maxYears)&&c.maxYears==fix(c.maxYears)&&c.maxYears>=1&&c.maxYears<=60);
assert(c.consecutiveYears>=2&&c.consecutiveYears==fix(c.consecutiveYears));
assert(c.productionYears==1,'One confirmation/production year is required.');
assert(all([c.tolT c.tolH c.tolWater c.tolIce c.roiTolT c.roiTolH ...
 c.roiTolWater c.roiTolIce c.roiRadius c.tolLeakRelative c.tolLeakAbsolute]>0));
if ~isfolder(c.outputFolder), mkdir(c.outputFolder); end
checkpointFile=fullfile(c.outputFolder,'spinup_progress.mat');
activeCheckpointFile=fullfile(c.outputFolder,'active_year_checkpoint.mat');
o=struct('meshFile',fullfile(root,'example_mesh_repaired_interface.mat'),'durationDays',365, ...
 'outputEveryDays',365/12,'plot',false,'verbose',false,'saveOutput',true);
o.initial=struct('mode','profile','surfaceT',-1,'gradient',.006);
o.mechanics=struct('enable',false); % no feedback to TH; skip repeated elastic solves
o.leakage=struct('segmentCount',6,'segmentStartAngleRad',0, ...
 'interiorHead',0,'jointAnglesRad',[], ...
 'jointWidth_m',[],'jointConductance',[],'jointEntryHead',[], ...
 'jointFreeze',joint_freeze_defaults());
o.dtInitial=.1*86400; o.dtMax=86400; o.dtMin=86400/16384;
% tol=2e-8 remains far below the physical stopping tolerances and avoids
% demanding a residual below the finite-difference Jacobian noise floor.
o.solver=struct('fd',1e-8,'maxIter',60,'tol',2e-8, ...
 'jacobianReuse',4,'jacobianRefreshRatio',.65, ...
 'dtTargetNewtonFraction',.4,'dtAcceptedGrowthFloor',1.05, ...
 'dtHighIterationFraction',.75,'dtGrowthSafety',1.10, ...
 'dtGrowthMax',1.8);
o=combine(o,c.model);
assert(~isfield(o,'yearDays')||o.yearDays==365);
assert(~isfield(o,'day')||o.day==86400);
assert(o.dtMax<=4*86400&&o.dtMax>0, ...
 'Maximum step must not exceed four days; the joint-freezing case needs its own convergence study.');
assert(any(strcmp(o.initial.mode,{'profile','fields'})), ...
 'Fresh runs use profile or measured fields; use resume for continuation.');
o.durationDays=365; o.outputEveryDays=365/12; o.plot=false;
o.outputPrefix=fullfile(c.outputFolder,'active_year');
columns={'year','rmsDeltaT_C','rmsDeltaH_m','rmsDeltaWater','rmsDeltaIce', ...
 'roiMaxDeltaT_C','roiMaxDeltaH_m','roiMaxDeltaWater','roiMaxDeltaIce', ...
 'globalMaxDeltaT_C','globalMaxDeltaH_m','globalMaxDeltaWater','globalMaxDeltaIce', ...
 'leak_m3_per_m','deltaAnnualLeak_m3_per_m','relativeAnnualLeak', ...
 'relativeDailyLeakL1','maxSegmentAnnualRelative','maxSegmentDailyRelative', ...
 'maxJointAnnualRelative','maxJointDailyRelative', ...
 'waterDefect_m3_per_m','energyDefect_J_per_m', ...
 'sumAbsWaterDefect','sumAbsEnergyDefect','acceptedSteps','maxResidual','pass'};
previous=[]; diagnostics=[]; legacyDiagnostics=[]; consecutive=0; convergedYear=NaN; completed=0;
migrationNote='none'; restartNote='none'; numericalFailure='';
if c.resume && ~isfile(checkpointFile) && isfile(activeCheckpointFile)
 % There is no completed-year state with which to assemble year 1.
 c.resume=false;
end
if c.resume
 assert(isfile(checkpointFile),'No completed-year checkpoint exists.');
 d=load(checkpointFile,'progress'); z=d.progress;
 previous=z.lastYear; completed=z.completedYears;
 assert(isfield(z,'schemaVersion') && z.schemaVersion==7, ...
  'This v12 model requires a v12 annual checkpoint. Start a new output folder.');
 assert(strcmp(previous.version,'TH_v12_20260930_joint_freeze'), ...
  'This annual checkpoint is not from the v12 solver. Start a new output folder.');
 if z.schemaVersion==7
  expected=rmfield(c,{'resume','outputFolder'}); actual=rmfield(z.config,{'resume','outputFolder'});
  assert(isequaln(expected,actual),'Resume must use identical model, criteria and year limit.');
  diagnostics=z.diagnostics; consecutive=z.consecutive; convergedYear=z.convergedYear;
  if isfield(z,'legacyDiagnostics'), legacyDiagnostics=z.legacyDiagnostics; end
  if isfield(z,'restartNote'), restartNote=z.restartNote; end
 end
else
 assert(~isfile(checkpointFile),'Output already exists. Use resume=true or a new outputFolder.');
end
if completed==0 && isfile(activeCheckpointFile)
 restartNote='Year-1 monthly checkpoint ignored: no completed prior year; reran year 1 from initial fields.';
 fprintf('%s\n',restartNote);
 activeResume=struct('use',false);
else
 activeResume=loadActiveYear(c.outputFolder,completed,o.leakage);
 if isfile(activeCheckpointFile) && ~activeResume.use && ...
         ~strcmp(activeResume.reason,'no compatible active-year checkpoint')
  fprintf('Ignoring active-year checkpoint: %s. Year %d will restart.\n', ...
   activeResume.reason,completed+1);
 end
end
if activeResume.use
 fprintf(['Found valid v12 month checkpoint at day %.8g. ' ...
  'Year %d will continue from this state.\n'],activeResume.day,completed+1);
end
status='maximum_years_reached_unconverged';
lastYearToRun=c.maxYears;
if ~isnan(convergedYear)&&completed>=convergedYear+c.productionYears
 lastYearToRun=completed;
end
for year=completed+1:lastYearToRun
 usingPartial=activeResume.use && year==completed+1;
 if usingPartial
  assert(~isempty(previous),'An active-year checkpoint requires a completed previous year.');
  o.mesh=previous.mesh;
  o.initial=struct('mode','checkpoint','checkpointFile',activeResume.file);
  o.durationDays=year*365-activeResume.day;
  savePartialPrefix(c.outputFolder,year,activeResume);
 elseif ~isempty(previous)
  o.mesh=previous.mesh; o.initial=struct('mode','fields','T',previous.T, ...
   'h',previous.h,'startDay',(year-1)*365);
  o.durationDays=365;
 else
  o.initial.startDay=0;
  o.durationDays=365;
 end
 fprintf('Starting year %d/%d (total physical-time cap includes production).\n',year,c.maxYears);
 try
  out=TH_2D_segment_v12(o);
 catch ME
  if ~strcmp(ME.identifier,'TH:nonconvergence'), rethrow(ME); end
  numericalFailure=sprintf('Year %d stopped at the last accepted substep: %s',year,ME.message);
  warning('TH:annualFailure','%s',numericalFailure); status='numerical_failure_saved'; break
 end
 if usingPartial
  assert(size(activeResume.hist,2)==size(out.history,2), ...
   'Partial-year history schema differs from v12 history.');
  out.history=[activeResume.hist;out.history];
  out.leakage.segmentRate_m2_per_s=[activeResume.segmentHist;out.leakage.segmentRate_m2_per_s];
  out.leakage.jointRate_m2_per_s=[activeResume.jointHist;out.leakage.jointRate_m2_per_s];
  out.leakage.bodyRate_m2_per_s=[activeResume.bodyHist;out.leakage.bodyRate_m2_per_s];
  out.leakage.jointFreeze=freezeHistory( ...
   [activeResume.jointFreezeHist;out.leakage.jointFreeze.history],out.parameters.leakage);
  out.leakage=summarize_segment_leakage(out.history,out.leakage,(year-1)*365);
  out.snapshots=[activeResume.snaps,out.snapshots];
  out.initialX=[previous.T;previous.h];
  out.partialYearResume=struct('used',true,'resumeDay',activeResume.day, ...
   'phaseSnapshotsComplete',numel(out.snapshots)==12);
  out.solverWork=workFromHistory(out.history,out.solverWork.rejectedSteps);
  save([o.outputPrefix '_results.mat'],'out','-v7');
 else
  out.partialYearResume=struct('used',false,'resumeDay',NaN, ...
   'phaseSnapshotsComplete',numel(out.snapshots)==12);
 end
 assert(abs(out.timeDay-year*365)<1e-7&&all(isfinite(out.history(:))));
 assert(all(diff(out.history(:,1))>0)&&abs(sum(out.history(:,2))-365)<1e-6, ...
  'Annual history is incomplete or contains duplicated time steps.');
 assert(validFreezeHistory(out.leakage.jointFreeze.history,size(out.history,1), ...
  out.parameters.leakage),'Incomplete or invalid annual joint-freezing history.');
 h=out.history; daily=out.leakage.dailyTotalVolume_m3_per_m;
 volume=out.leakage.totalVolume_m3_per_m;
 metric=inf(1,12); dv=inf; rv=inf; rq=inf; dailyDifference=inf; oldVolume=0;
 segmentAnnualRelative=inf; segmentDailyRelative=inf;
 jointAnnualRelative=inf; jointDailyRelative=inf;
 segmentPass=false; jointPass=false;
 if ~isempty(previous)
  if numel(out.snapshots)==12 && numel(previous.snapshots)==12
   metric=phaseDifference(out,previous,c);
  else
   warning('TH:incompletePhaseHistory', ...
    'Year %d has incomplete monthly phase history; field criteria fail safely for this comparison.',year);
  end
  oldDaily=previous.leakage.dailyTotalVolume_m3_per_m;
  oldVolume=previous.leakage.totalVolume_m3_per_m;
  dv=abs(volume-oldVolume); denominator=max([volume oldVolume c.tolLeakAbsolute]);
  dailyDifference=sum(abs(daily-oldDaily));
  rv=dv/denominator; rq=dailyDifference/denominator;
  [segmentPass,segmentAnnualRelative,segmentDailyRelative]=leakGroupPass( ...
   out.leakage.dailyVolumeBySegment_m3_per_m, ...
   previous.leakage.dailyVolumeBySegment_m3_per_m,c);
  [jointPass,jointAnnualRelative,jointDailyRelative]=leakGroupPass( ...
   out.leakage.dailyVolumeByJoint_m3_per_m, ...
   previous.leakage.dailyVolumeByJoint_m3_per_m,c);
 end
 pass=all(metric(1:4)<=[c.tolT c.tolH c.tolWater c.tolIce]) ...
  && all(metric(5:8)<=[c.roiTolT c.roiTolH c.roiTolWater c.roiTolIce]) ...
  && dv<=c.tolLeakAbsolute+c.tolLeakRelative*max(volume,oldVolume) ...
  && dailyDifference<=c.tolLeakAbsolute+c.tolLeakRelative*max(volume,oldVolume) ...
  && segmentPass && jointPass;
 if pass, consecutive=consecutive+1; else, consecutive=0; convergedYear=NaN; end
 if consecutive>=c.consecutiveYears&&isnan(convergedYear), convergedYear=year; end
 row=[year metric volume dv rv rq segmentAnnualRelative segmentDailyRelative ...
  jointAnnualRelative jointDailyRelative sum(h(:,8)) sum(h(:,9)) ...
  sum(abs(h(:,8))) sum(abs(h(:,9))) size(h,1) max(h(:,11)) double(pass)];
 diagnostics=[diagnostics;row]; %#ok<AGROW>
 fprintf(['Year %d: RMS[dT dh dW dI]=[%.3g %.3g %.3g %.3g], ' ...
  'ROImax=[%.3g %.3g %.3g %.3g]; globalMax=[%.3g %.3g %.3g %.3g]; ' ...
  'leak=%.6g L/m; annual=%.3g%%; dailyL1=%.3g%%; consecutive=%d\n'], ...
  year,metric,1000*volume,100*rv,100*rq,consecutive);
 previous=out; completed=year;
 activeResume.use=false;
 progress=struct('schemaVersion',7,'config',c,'completedYears',completed,'lastYear',out, ...
  'diagnostics',diagnostics,'columns',{columns},'consecutive',consecutive,'convergedYear',convergedYear);
 progress.legacyDiagnostics=legacyDiagnostics; progress.migrationNote=migrationNote;
 progress.restartNote=restartNote;
 % Atomic completed-year checkpoint; an interruption loses at most the active year.
 tempFile=fullfile(c.outputFolder,'spinup_progress_tmp.mat');
 save(tempFile,'progress','-v7'); movefile(tempFile,checkpointFile,'f');
 if ~isnan(convergedYear)&&year>=convergedYear+c.productionYears
  status='converged_and_confirmation_year_completed'; break
 end
 % Optional user-requested orderly pause; remove this marker before resuming.
 if isfile(fullfile(c.outputFolder,'PAUSE_AFTER_YEAR'))
  status='paused_at_year_boundary'; break
 end
end
if ~any(strcmp(status,{'converged_and_confirmation_year_completed','paused_at_year_boundary', ...
 'numerical_failure_saved'}))&&~isnan(convergedYear)
 if completed>=convergedYear+1
  status='converged_and_confirmation_year_completed';
 else
  status='converged_at_cap_no_confirmation_year';
 end
end
result=struct('status',status,'completedYears',completed,'convergedYear',convergedYear, ...
 'periodicConverged',~isnan(convergedYear),'productionComplete',strcmp(status,'converged_and_confirmation_year_completed'), ...
 'diagnostics',diagnostics,'columns',{columns},'lastYear',previous,'config',c);
result.schemaVersion=7; result.legacyDiagnostics=legacyDiagnostics;
result.migrationNote=migrationNote; result.numericalFailure=numericalFailure;
result.restartNote=restartNote;
save(fullfile(c.outputFolder,'TH_segment_spinup_results.mat'),'result','-v7');
fprintf('STOP: %s; %d years; never exceeds %d total years.\n',status,completed,c.maxYears);
end

function r=loadActiveYear(folder,completed,leakConfig)
% Recover a v12 monthly checkpoint newer than the last complete year.
emptySnaps=struct('timeDay',{},'T',{},'h',{},'thetaW',{},'thetaI',{});
r=struct('use',false,'file','','day',NaN,'hist',zeros(0,14), ...
 'snaps',{emptySnaps},'segmentHist',zeros(0,leakConfig.segmentCount), ...
 'jointHist',zeros(0,numel(leakConfig.jointAnglesRad)), ...
 'jointFreezeHist',zeros(0,4*numel(leakConfig.jointAnglesRad)), ...
 'bodyHist',zeros(0,1),'reason','no compatible active-year checkpoint');
file=fullfile(folder,'active_year_checkpoint.mat');
if ~isfile(file), return; end
try
z=load(file);
if ~all(isfield(z,{'checkpoint','hist','leakHistSegment','leakHistJoint', ...
        'leakHistBody','jointFreezeHist'}))
 r.reason='incomplete active-year checkpoint'; return
end
cp=z.checkpoint;
if ~isfield(cp,'version')||~strcmp(cp.version,'TH_v12_20260930_joint_freeze'), return; end
if ~isequaln(cp.parameters.leakage,leakConfig), return; end
day=cp.time/86400; year=completed+1; first=(year-1)*365; last=year*365;
if ~(day>first+1e-7 && day<last-1e-7), return; end
phase=(day-first)/(365/12);
if abs(phase-round(phase))>1e-6, return; end
if isempty(z.hist)||size(z.hist,2)~=14||abs(z.hist(end,1)-day)>1e-6, return; end
if size(z.leakHistSegment,1)~=size(z.hist,1) || ...
        size(z.leakHistSegment,2)~=leakConfig.segmentCount || ...
        size(z.leakHistJoint,1)~=size(z.hist,1) || ...
        size(z.leakHistJoint,2)~=numel(leakConfig.jointAnglesRad) || ...
        size(z.leakHistBody,1)~=size(z.hist,1) || ...
        ~validFreezeHistory(z.jointFreezeHist,size(z.hist,1),leakConfig)
    return;
end
hist=z.hist; snaps=emptySnaps;
segmentHist=z.leakHistSegment; jointHist=z.leakHistJoint;
bodyHist=z.leakHistBody;
jointFreezeHist=z.jointFreezeHist;
if isfield(z,'snaps'), snaps=reshape(z.snaps,1,[]); end

% Repeated interruptions remain resumable: the prefix stores all history
% before the current engine segment, because the engine checkpoint itself
% contains history only since its most recent start.
prefixFile=fullfile(folder,'partial_year_prefix.mat');
if isfile(prefixFile)
 q=load(prefixFile,'partialYearBase');
 if isfield(q,'partialYearBase')
  b=q.partialYearBase;
  if isfield(b,'year')&&b.year==year&&isfield(b,'checkpointDay')
   assert(isfield(b,'jointFreezeHist'), ...
    'Partial-year prefix lacks v12 joint-freezing history.');
   if abs(b.checkpointDay-day)<1e-6
    hist=b.hist; snaps=reshape(b.snaps,1,[]);
    segmentHist=b.segmentHist; jointHist=b.jointHist; bodyHist=b.bodyHist;
    jointFreezeHist=b.jointFreezeHist;
   elseif b.checkpointDay<day-1e-7 && z.hist(1,1)>b.checkpointDay-1e-6
    hist=[b.hist;z.hist]; snaps=[reshape(b.snaps,1,[]),snaps];
    segmentHist=[b.segmentHist;z.leakHistSegment];
    jointHist=[b.jointHist;z.leakHistJoint];
    bodyHist=[b.bodyHist;z.leakHistBody];
    jointFreezeHist=[b.jointFreezeHist;z.jointFreezeHist];
   end
  end
 end
end
assert(all(diff(hist(:,1))>0)&&abs(hist(end,1)-day)<1e-6, ...
 'Corrupt or duplicated partial-year time history.');
assert(size(segmentHist,1)==size(hist,1)&&size(jointHist,1)==size(hist,1) ...
 && size(bodyHist,1)==size(hist,1),'Incomplete partial-year leakage history.');
assert(validFreezeHistory(jointFreezeHist,size(hist,1),leakConfig), ...
 'Incomplete or invalid partial-year joint-freezing history.');
r=struct('use',true,'file',file,'day',day,'hist',hist,'snaps',{snaps}, ...
 'segmentHist',segmentHist,'jointHist',jointHist,'bodyHist',bodyHist, ...
 'jointFreezeHist',jointFreezeHist, ...
 'reason','compatible month checkpoint');
catch ME
 r.use=false;
 r.reason=['unreadable or incomplete active-year checkpoint: ' ME.message];
end
end

function savePartialPrefix(folder,year,r)
partialYearBase=struct('year',year,'checkpointDay',r.day, ...
 'hist',r.hist,'snaps',{r.snaps},'segmentHist',r.segmentHist, ...
 'jointHist',r.jointHist,'bodyHist',r.bodyHist, ...
 'jointFreezeHist',r.jointFreezeHist); %#ok<NASGU>
tmp=fullfile(folder,'partial_year_prefix_tmp.mat');
file=fullfile(folder,'partial_year_prefix.mat');
save(tmp,'partialYearBase','-v7'); movefile(tmp,file,'f');
end

function valid=validFreezeHistory(history,numberSteps,leakConfig)
nj=numel(leakConfig.jointAnglesRad);
valid=isequal(size(history),[numberSteps,4*nj]);
if ~valid, return; end
factor=history(:,1:nj);
minimum=history(:,nj+(1:nj));
maximum=history(:,2*nj+(1:nj));
frozen=history(:,3*nj+(1:nj));
valid=all(isfinite(factor(:)))&&all(factor(:)>=0) ...
 && all(isfinite(frozen(:)))&&all(frozen(:)>=0)&&all(frozen(:)<=1);
if leakConfig.jointFreeze.enabled
 valid=valid&&all(isfinite(minimum(:)))&&all(isfinite(maximum(:))) ...
  && all(minimum(:)<=maximum(:));
else
 % Disabled records deliberately omit temperatures; NaN is not corruption.
 valid=valid&&all(isnan(minimum(:)))&&all(isnan(maximum(:))) ...
  && all(factor(:)==1)&&all(frozen(:)==0);
end
end

function f=freezeHistory(history,leakConfig)
nj=numel(leakConfig.jointAnglesRad);
c=leakConfig.jointConductance(:)';
if nj==0, c=zeros(1,0); elseif isscalar(c), c=repmat(c,1,nj); end
f=struct('enabled',leakConfig.jointFreeze.enabled,'history',history, ...
 'factor',history(:,1:nj),'minTemperature_C',history(:,nj+(1:nj)), ...
 'maxTemperature_C',history(:,2*nj+(1:nj)), ...
 'frozenLengthFraction',history(:,3*nj+(1:nj)), ...
 'conductance_m_per_s',history(:,1:nj).*c);
end

function w=workFromHistory(h,rejected)
w=struct('acceptedSteps',size(h,1),'nonlinearIterations',sum(h(:,3)), ...
 'jacobianBuilds',sum(h(:,12)),'linearSolves',sum(h(:,13)), ...
 'residualEvaluations',sum(h(:,14)),'rejectedSteps',rejected);
end

function [passed,annualRelative,dailyRelative]=leakGroupPass(daily,currentPrevious,c)
assert(isequal(size(daily),size(currentPrevious)));
n=size(daily,2);
if n==0
 passed=true; annualRelative=0; dailyRelative=0; return
end
current=sum(daily,1); previous=sum(currentPrevious,1);
annualDelta=abs(current-previous);
dailyL1=sum(abs(daily-currentPrevious),1);
absolute=c.tolLeakAbsolute/n;
limit=absolute+c.tolLeakRelative*max(current,previous);
denominator=max(max(current,previous),absolute);
passed=all(annualDelta<=limit) && all(dailyL1<=limit);
annualRelative=max(annualDelta./denominator);
dailyRelative=max(dailyL1./denominator);
end

function a=combine(a,b)
f=fieldnames(b);
for k=1:numel(f)
    if isfield(a,f{k})&&isstruct(a.(f{k}))&&isstruct(b.(f{k}))
        a.(f{k})=combine(a.(f{k}),b.(f{k}));
    else
        a.(f{k})=b.(f{k});
    end
end
end

function m=phaseDifference(a,b,c)
% Compare 12 matching phases AND year endpoints across all nodes.
assert(numel(a.snapshots)==12&&numel(b.snapshots)==12);
m=zeros(1,12);
for k=1:12
 x=a.snapshots(k); y=b.snapshots(k);
 assert(abs((x.timeDay-y.timeDay)-365)<1e-6);
 m=max(m,oneDifference(x.T,x.h,x.thetaW+x.thetaI,x.thetaI, ...
  y.T,y.h,y.thetaW+y.thetaI,y.thetaI,a,c));
end
m=max(m,oneDifference(a.T,a.h,a.state.total,a.state.thetaI, ...
 b.T,b.h,b.state.total,b.state.thetaI,a,c));
end

function m=oneDifference(T,h,w,i,T0,h0,w0,i0,out,c)
d={abs(T-T0),abs(h-h0),abs(w-w0),abs(i-i0)};
rmsv=cellfun(@(x)sqrt(mean(x.^2)),d);
xy=out.pAll;
if isfield(out.mesh,'xc')&&isfield(out.mesh,'yc')
 center=[out.mesh.xc out.mesh.yc];
else
 center=(min(xy,[],1)+max(xy,[],1))/2;
end
roiT=hypot(xy(:,1)-center(1),xy(:,2)-center(2))<=c.roiRadius;
soil=out.soilNodes; roiH=roiT(soil);
roiv=[max(d{1}(roiT)),max(d{2}(roiH)),max(d{3}(roiH)),max(d{4}(roiH))];
globalv=cellfun(@max,d);
m=[rmsv roiv globalv];
end

