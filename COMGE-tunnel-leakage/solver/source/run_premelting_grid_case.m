function result = run_premelting_grid_case(folder,nAngular,nRadial,nJointCells,fine)
%RUN_PREMELTING_GRID_CASE Single immutable nominal verification trajectory.
% Run from a fresh MATLAB process or clear functions before this entry.
cfg=premelting_research_case(); cfg.caseID=sprintf('A%d_R%d_J%d_fine%d',nAngular,nRadial,nJointCells,fine);
cfg.nAngular=nAngular; cfg.nRadial=nRadial; cfg.nJointCells=nJointCells;
% Preserve the original constitutive baseline for the recorded mesh study.
% A viscosity-law comparison is a different physical case, not a refinement.
cfg.parameters.viscosityModel='constant';
cfg.runOptions.stepController='growth';
cfg.useCyclicSymmetry=true; cfg.runOptions.verbose=true;
if fine
    cfg.runOptions.dtDays=cfg.runOptions.dtDays/2;
    cfg.runOptions.maximumPhaseIncrement=cfg.runOptions.maximumPhaseIncrement/2;
    cfg.runOptions.maximumPressureIncrement_Pa=cfg.runOptions.maximumPressureIncrement_Pa/2;
    cfg.runOptions.maximumActiveIcePressureIncrement_Pa=cfg.runOptions.maximumActiveIcePressureIncrement_Pa/2;
end
result=run_premelting_sample(cfg,folder);
fprintf('Sample %s status: %s\n',cfg.caseID,result.status);
if strcmp(result.status,'solver_completed'), disp(result.qoi);
else, error('PhaseSwitch:CampaignFailed','%s: %s',result.failure.identifier,result.failure.message); end
end
