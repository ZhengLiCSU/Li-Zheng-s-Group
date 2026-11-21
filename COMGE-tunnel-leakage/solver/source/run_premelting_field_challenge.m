function result = run_premelting_field_challenge(folder,nCells)
%RUN_PREMELTING_FIELD_CHALLENGE Full eight-joint deterministic nonuniform case.
% No distribution or reliability estimate is assigned. A 3-cell run checks
% the material-field gateway; it does NOT establish spatial convergence.
arguments
    folder (1,:) char
    nCells (1,1) double {mustBeInteger,mustBePositive} = 3
end
cfg=premelting_research_case(); cfg.caseID=sprintf('nonuniform_full8_J%d',nCells);
cfg.nAngular=16; cfg.nRadial=nCells; cfg.nJointCells=nCells;
cfg.useCyclicSymmetry=false; cfg.runOptions.verbose=true;
xi=linspace(0,1,5); joint=(0:7)';
values=5e-6*exp(.12*cos(2*pi*joint/8)+.1*sin(pi*xi));
cfg.apertureField=struct('xi',xi,'aperture_m',values);
result=run_premelting_sample(cfg,folder);
assert(strcmp(result.status,'solver_completed'),'PhaseSwitch:CampaignFailed', ...
    'The challenge failed; inspect the retained sample status and trajectory.');
disp(result.qoi);
end
