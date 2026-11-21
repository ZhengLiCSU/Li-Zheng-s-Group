function checks=verify_segment_v11()
%VERIFY_SEGMENT_V11 Check total-head leakage, budgets, and controller resume.
checks.physics=TH_2D_segment_v11('selftest');
rate=2e-8; h=zeros(365,14); h(:,1)=(1:365)'; h(:,2)=1; h(:,10)=rate;
leak=struct('segmentRate_m2_per_s',repmat([.3 .7]*rate,365,1), ...
    'jointRate_m2_per_s',zeros(365,0), ...
    'bodyRate_m2_per_s',rate*ones(365,1));
leak=summarize_segment_leakage(h,leak,0);
assert(abs(leak.totalVolume_m3_per_m-rate*365*86400)<1e-12);
assert(abs(sum(leak.monthlyTotalVolume_m3_per_m)-leak.totalVolume_m3_per_m)<1e-12);
assert(size(leak.dailyVolumeBySegment_m3_per_m,1)==365);
base=tempname; mkdir(base); clean=onCleanup(@() rmdir(base,'s')); %#ok<NASGU>
S=struct('pAll',[0 0;1 0;0 1;1 1],'tSoil',[1 2 3;2 4 3]);
m=struct('mesh',S,'topMin',2,'topMax',2,'topHTC',inf,'geothermalFlux',0, ...
    'hydroFar',true,'waterTableDepth',-1,'topSeepage',false,'tunnelDrain',false, ...
    'saveOutput',false,'initial',struct('mode','profile','surfaceT',2,'gradient',0));
forced=m;
forced.durationDays=1; forced.outputEveryDays=1;
forced.verbose=false; forced.plot=false;
forced.forcing=struct('day',[0 1 365], ...
    'topAirC',[2 4 2],'waterTableDepth_m',[-1 -2 -1]);
f=TH_2D_segment_v11(forced);
assert(max(abs(f.T(S.pAll(:,2)==0)-4))<1e-10);
assert(abs(f.h(4)-(S.pAll(4,2)+2))<1e-9);
realMesh=TH_2D_segment_v11(struct('durationDays',1,'outputEveryDays',1, ...
    'saveOutput',false,'plot',false,'verbose',false, ...
    'mechanics',struct('enable',false)));
assert(realMesh.numerics.negativeLiningEdges==0 && ...
    realMesh.numerics.negativeHydraulicEdges==8, ...
    'The v11 default mesh must have no negative lining hydraulic edges.');
c=struct('model',m,'maxYears',6,'outputFolder',fullfile(base,'stationary'));
a=run_TH_segment_spinup_v11(c);
assert(a.convergedYear==4&&a.completedYears==5&&a.productionComplete);
assert(size(a.lastYear.leakage.monthlyVolumeBySegment_m3_per_m,1)==12);
c.resume=true; b=run_TH_segment_spinup_v11(c);
assert(isequaln(a.diagnostics,b.diagnostics));
% Resume the second year from a real month checkpoint and preserve its rates.
c=struct('model',m,'maxYears',2,'outputFolder',fullfile(base,'monthly'));
mkdir(c.outputFolder);
marker=fullfile(c.outputFolder,'PAUSE_AFTER_YEAR');
fid=fopen(marker,'w'); fclose(fid);
paused=run_TH_segment_spinup_v11(c);
assert(paused.completedYears==1);
delete(marker);
opt=paused.lastYear.parameters;
opt.initial.mode='fields'; opt.initial.T=paused.lastYear.T;
opt.initial.h=paused.lastYear.h; opt.initial.startDay=365;
opt.durationDays=365/12; opt.outputEveryDays=365/12;
opt.outputPrefix=fullfile(c.outputFolder,'active_year');
opt.saveOutput=true;
partial=TH_2D_segment_v11(opt);
assert(abs(partial.timeDay-(365+365/12))<1e-8);
% A checkpoint can contain retries from earlier in the same year. Their
% count must survive the annual controller's monthly restart.
checkpointFile=fullfile(c.outputFolder,'active_year_checkpoint.mat');
checkpointData=load(checkpointFile,'checkpoint');
checkpoint=checkpointData.checkpoint;
checkpoint.rejectedSteps=7;
save(checkpointFile,'checkpoint','-append');
c.resume=true; resumed=run_TH_segment_spinup_v11(c);
assert(resumed.completedYears==2&&resumed.lastYear.partialYearResume.used);
assert(size(resumed.lastYear.history,1)==size(resumed.lastYear.leakage.segmentRate_m2_per_s,1));
assert(size(resumed.lastYear.leakage.dailyTotalVolume_m3_per_m,1)==365);
assert(resumed.lastYear.solverWork.rejectedSteps>=7, ...
    'The annual restart discarded the monthly checkpoint retry count.');
assert(resumed.lastYear.solverWork.acceptedSteps==size(resumed.lastYear.history,1));
% With no completed prior year, a first-year monthly checkpoint cannot be
% assembled into a full annual record and must trigger a clean rerun.
c=struct('model',m,'maxYears',1,'outputFolder',fullfile(base,'first_year_restart'), ...
    'resume',true);
mkdir(c.outputFolder);
fid=fopen(fullfile(c.outputFolder,'active_year_checkpoint.mat'),'w');
assert(fid>0); fwrite(fid,'interrupted first year'); fclose(fid);
firstYear=run_TH_segment_spinup_v11(c);
assert(firstYear.completedYears==1 && contains(firstYear.restartNote, ...
    'Year-1 monthly checkpoint ignored'));
% An unreadable active checkpoint after a completed year is also ignored.
c=struct('model',m,'maxYears',2,'outputFolder',fullfile(base,'corrupt_month_restart'));
mkdir(c.outputFolder);
marker=fullfile(c.outputFolder,'PAUSE_AFTER_YEAR');
fid=fopen(marker,'w'); assert(fid>0); fclose(fid);
paused=run_TH_segment_spinup_v11(c);
assert(paused.completedYears==1); delete(marker);
fid=fopen(fullfile(c.outputFolder,'active_year_checkpoint.mat'),'w');
assert(fid>0); fwrite(fid,'interrupted monthly save'); fclose(fid);
c.resume=true;
restarted=run_TH_segment_spinup_v11(c);
assert(restarted.completedYears==2 && ~restarted.lastYear.partialYearResume.used);
checks.annualIntegration=true;
checks.observedForcing=true;
checks.repairedMesh=true;
checks.stationaryAndResume=true;
checks.monthlyCheckpointResume=true;
checks.rejectedStepResume=true;
checks.firstYearRestart=true;
checks.corruptMonthRestart=true;
disp(checks);
end

