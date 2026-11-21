function result = run_premelting_sample(cfg,folder)
%RUN_PREMELTING_SAMPLE Deterministic realization gateway for later UQ.
% Failed inputs/solves are retained with no substituted zero leakage.
arguments
    cfg (1,1) struct
    folder (1,:) char
end
assert(~isfolder(folder) && ~isfile(folder),'PhaseSwitch:OutputExists','Use a fresh sample folder; evidence is never overwritten.');
sampleTimer=tic;
mkdir(folder); result=struct('status','pending','qoi',[],'failure',[]);
save(fullfile(folder,'input.mat'),'cfg','-v7');
sourceFolder=fullfile(folder,'source'); mkdir(sourceFolder);
base=fileparts(mfilename('fullpath')); files=dir(fullfile(base,'*.m'));
for k=1:numel(files), copyfile(fullfile(base,files(k).name),sourceFolder); end
manifest=premelting_source_manifest(sourceFolder);
save(fullfile(folder,'source_manifest.mat'),'manifest','-v7');
environment=struct('MATLABVersion',version,'release',version('-release'), ...
    'computer',computer,'computationThreads',maxNumCompThreads);
save(fullfile(folder,'environment.mat'),'environment','-v7');
try
    expected=fieldnames(premelting_research_case());
    assert(isequal(sort(fieldnames(cfg)),sort(expected)) && cfg.schema==1,'PhaseSwitch:Input','Unsupported or incomplete input schema.');
    assert(ischar(cfg.caseID) && ~isempty(cfg.caseID),'PhaseSwitch:Input','A case ID is required.');
    assert(~cfg.cohesiveParameters.enabled,'PhaseSwitch:Qualification','Damage-enabled switching is not qualified by this component gateway.');
    assert(~isfield(cfg.runOptions,'saveFile'),'PhaseSwitch:Input','The gateway owns the evidence path.');
    m=build_cohesive_threephase_ring(make_refined_thm_verification_mesh(cfg.nAngular,cfg.nRadial),cfg.nJointCells,cfg.parameters,cfg.cohesiveParameters);
    meshEvidence=check_premelting_field_mesh(m);
    p=cfg.parameters; fieldMetadata=[];
    if ~isempty(cfg.apertureField), [p,fieldMetadata]=apply_premelting_aperture_field(m,p,cfg.apertureField); end
    if cfg.useCyclicSymmetry, [m,p]=reduce_premelting_cyclic_symmetry(m,p); end
    save(fullfile(folder,'assembled_input.mat'),'m','p','fieldMetadata','meshEvidence','-v7');
    options=cfg.runOptions; options.saveFile=fullfile(folder,'trajectory.mat');
    out=solve_premelting_switch_ring(m,p,cfg.filmParameters,options);
    result.status='solver_completed'; result.qoi=summarize_premelting_sample(out);
catch exception
    result.status='numerical_failure';
    if contains(exception.identifier,{'Domain','Field','Input','Shape','Qualification','Symmetry','Material','Geometry','TimeControl','TraceResolution'})
        result.status='input_or_physical_domain_failure';
    end
    result.failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
end
save(fullfile(folder,'sample_result.mat'),'result','-v7');
fid=fopen(fullfile(folder,'status.json'),'w'); assert(fid>=0,'PhaseSwitch:IO','Cannot save sample status.');
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(result));
clear cleanup
% Includes source archival, assembly, solve, QoI reconstruction and result
% commits. Only writing this small timing artifact itself is outside it.
runtime=struct('endToEndWallSeconds',toc(sampleTimer), ...
    'scope','Gateway invocation through result commit; MATLAB startup excluded.');
save(fullfile(folder,'runtime.mat'),'runtime','-v7');
end
