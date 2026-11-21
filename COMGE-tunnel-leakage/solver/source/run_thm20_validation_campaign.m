function result=run_thm20_validation_campaign(resume)
%RUN_THM20_VALIDATION_CAMPAIGN Full default-mesh 20-year THM verification.
% No cycle shortcut, imported TH state, changed physical parameter, or
% reduced phase quadrature is used. Numerical checks precede a fresh run.
if nargin<1, resume=false; end
root=fileparts(mfilename('fullpath')); addpath(root); maxNumCompThreads(4);
folder=fullfile(root,'thm20_20260930_analytic_p1');
if ~isfolder(folder), mkdir(folder); end
if ~resume
    tests=runtests({'TestTHMBulkConstitutive','TestJointOpeningP1Hydraulics', ...
        'TestV16THM','TestTHMAnnualContinuation'});
    assertSuccess(tests);
    save(fullfile(folder,'preflight_tests.mat'),'tests');
    sourceFolder=fullfile(folder,'source_at_start'); mkdir(sourceFolder);
    sources=dir(fullfile(root,'*.m'));
    for k=1:numel(sources), copyfile(fullfile(root,sources(k).name),sourceFolder); end
    startedAt=datetime('now','TimeZone','UTC');
    save(fullfile(folder,'campaign_started.mat'),'startedAt');
end
cfg=struct('maxYears',20,'resume',logical(resume),'outputFolder',folder, ...
    'model',struct('solver',struct('jacobian','analytic','reuseAcrossSteps',true), ...
        'integration',struct('soilQuadratureOrder',5)));
result=run_TH_segment_year20_thm(cfg);
if result.completedYears==20
    validation=validate_thm_20year(folder,struct('saveReport',true));
    assert(validation.passed,'THM20:Validation','Full-history validation did not pass.');
    save(fullfile(folder,'completed_validation.mat'),'validation');
    timeStepValidation=verify_thm20_time_step(result);
    save(fullfile(folder,'completed_target_year_time_step.mat'),'timeStepValidation');
    assert(timeStepValidation.passed,'THM20:TimeStepValidation', ...
        'Year-20 half-day leakage verification did not pass; all results have been retained.');
end
end
