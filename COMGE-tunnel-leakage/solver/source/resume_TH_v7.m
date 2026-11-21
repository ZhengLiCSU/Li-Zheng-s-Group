function result=resume_TH_v7()
% One-command safe continuation from the latest complete year/month checkpoint.
root=fileparts(mfilename('fullpath')); addpath(root);
folder=fullfile(root,'spinup60_output');
assert(isfile(fullfile(folder,'spinup_progress.mat')), ...
    'spinup_progress.mat is missing; copy the existing spinup60_output folder here.');
logFile=fullfile(root,'runlog_v7_spinup.txt');
diary(logFile); finishLog=onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('\n=== v7 continuation started %s ===\n',datestr(now,31));
result=run_TH_spinup60(struct('resume',true));
save(fullfile(root,'TH_spinup60_v7_result.mat'),'result','-v7.3');
fprintf('=== v7 continuation finished %s: %s ===\n',datestr(now,31),result.status);
end
