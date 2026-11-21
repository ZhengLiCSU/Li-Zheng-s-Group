function result=run_joint_freeze_case(name)
%RUN_JOINT_FREEZE_CASE Reproducible full-history freeze mechanism contrasts.
root=fileparts(mfilename('fullpath')); addpath(root);
assert(any(strcmp(name,{'coldIce','coldViscosity','warmIce'})));
folder=fullfile(root,'joint_freeze_study_20260930',name);
model=struct('dtMax',4*86400);
if strcmp(name,'coldViscosity')
    model.leakage=struct('jointFreeze',struct('iceEnabled',false));
elseif strcmp(name,'warmIce')
    % Translate all temperatures by 10 C; retain the same hydraulic forcing.
    base=joint_only_th_model();
    model.topMin=base.topMin+10; model.topMax=base.topMax+10;
    model.tunnelMin=base.tunnelMin+10; model.tunnelMax=base.tunnelMax+10;
    model.groundwaterInflowT=base.groundwaterInflowT+10;
    model.initial=struct('surfaceT',base.initial.surfaceT+10);
end
maxNumCompThreads(4);
user=struct('model',model,'outputFolder',folder, ...
    'resume',isfile(fullfile(folder,'spinup_progress.mat')));
try
    result=run_TH_segment_year20_jointfreeze(user);
    assert(result.completedYears==20,'Full 20-year case did not finish.');
    fid=fopen(fullfile(folder,'COMPLETE.txt'),'w');
    fprintf(fid,'Q20_L_per_m=%.12g\nruntime_s=%.9g\n', ...
        result.year20AnnualLeakage_L_per_m,result.elapsedWallTime_s);
    fclose(fid);
catch ME
    if ~isfolder(folder), mkdir(folder); end
    fid=fopen(fullfile(folder,'ERROR.txt'),'w');
    fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    rethrow(ME);
end
end
