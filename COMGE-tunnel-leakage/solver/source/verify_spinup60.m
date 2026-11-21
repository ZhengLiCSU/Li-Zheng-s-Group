function checks=verify_spinup60()
% Exercise stopping and restart logic using an exact stationary TH state.
base=tempname; mkdir(base); clean=onCleanup(@() rmdir(base,'s')); %#ok<NASGU>
S=struct('pAll',[0 0;1 0;0 1;1 1],'tSoil',[1 2 3;2 4 3]);
m=struct('mesh',S,'topMin',2,'topMax',2,'topHTC',inf,'geothermalFlux',0, ...
 'hydroFar',true,'waterTableDepth',-1,'topSeepage',false,'tunnelDrain',false, ...
 'saveOutput',false,'initial',struct('mode','profile','surfaceT',2,'gradient',0));
c=struct('model',m,'maxYears',6,'outputFolder',fullfile(base,'stationary'));
a=run_TH_spinup60(c);
assert(a.convergedYear==4&&a.completedYears==5&&a.productionComplete);
c.resume=true; b=run_TH_spinup60(c);
assert(isequaln(a.diagnostics,b.diagnostics)&&b.completedYears==5);
c.resume=false; c.maxYears=2; c.outputFolder=fullfile(base,'cap2');
d=run_TH_spinup60(c);
assert(d.completedYears==2&&~d.periodicConverged&&~d.productionComplete);
c.maxYears=4; c.outputFolder=fullfile(base,'cap4'); e=run_TH_spinup60(c);
assert(e.completedYears==4&&e.periodicConverged&&~e.productionComplete);
rejected=false;
try, run_TH_spinup60(struct('maxYears',61)); catch, rejected=true; end
assert(rejected);
% An explicit year-boundary pause and portable completed-year checkpoint.
c.maxYears=6; c.outputFolder=fullfile(base,'pause'); mkdir(c.outputFolder);
fid=fopen(fullfile(c.outputFolder,'PAUSE_AFTER_YEAR'),'w'); fclose(fid);
f=run_TH_spinup60(c);
assert(f.completedYears==1&&strcmp(f.status,'paused_at_year_boundary'));
destination=fullfile(base,'moved'); mkdir(destination);
copyfile(fullfile(c.outputFolder,'spinup_progress.mat'),destination);
c.outputFolder=destination; c.resume=true; g=run_TH_spinup60(c);
assert(g.completedYears==5&&g.productionComplete&&isequaln(a.diagnostics,g.diagnostics));
c.tolT=.03; rejected=false;
try, run_TH_spinup60(c); catch, rejected=true; end
assert(rejected);
checks=struct('stationaryConvergenceAndConfirmation',true,'completedResumeNoExtraYear',true, ...
 'unconvergedCapStops',true,'convergedAtCapNoExtraYear',true,'rejectMoreThan60',true, ...
 'orderlyPauseAndPortableResume',true,'rejectChangedResumeCriteria',true);
disp(checks);
end
