function report = run_nonisothermal_ring_verification(folder)
%RUN_NONISOTHERMAL_RING_VERIFICATION Freeze complete source and audit runs.
% Failure of the saturated regime is a REQUIRED retained outcome, not a
% successful cold-cycle result. Zero thermal expansion is explicitly a
% limiting benchmark, never relabelled as calibrated tunnel behaviour.
arguments
    folder (1,:) char = fullfile(pwd,'nonisothermal_ring_20261003_verified')
end
assert(~isfolder(folder),'ThermoInterface:OutputExists','Use a new output directory.');
mkdir(folder); sourceFolder=fullfile(folder,'source_at_start'); mkdir(sourceFolder);
files={'joint_nonisothermal_face.m','thermo_interface_state.m','thermo_interface_residual.m', ...
    'build_thermo_interface_ring.m','thermo_interface_boundary.m','thermo_interface_step.m', ...
    'solve_nonisothermal_ring.m','diagnose_thermo_interface_domain.m','TestNonisothermalInterface.m', ...
    'run_nonisothermal_ring_verification.m','joint_thermodynamic_parameters.m','joint_helmholtz_state.m', ...
    'joint_adiabatic_equilibrium.m','joint_closed_equilibrium.m','build_segment_mechanical_mesh.m', ...
    'solve_segment_contact_fem.m','build_joint_gap_profiles.m','build_joint_thermal_paths.m', ...
    'joint_opening_p1_hydraulics.m','make_thm_verification_mesh.m'};
for k=1:numel(files), copyfile(which(files{k}),fullfile(sourceFolder,files{k})); end
diary(fullfile(folder,'verification.log')); finish=onCleanup(@()diary('off')); %#ok<NASGU>
tests=runtests('TestNonisothermalInterface'); assert(all([tests.Passed]),'ThermoInterface:Tests','Unit tests failed.');
issues=cell(size(files));
for k=1:numel(files), issues{k}=checkcode(files{k},'-id'); end
save(fullfile(folder,'tests_and_code_review.mat'),'tests','files','issues','-v7');
p=joint_thermodynamic_parameters(); S=make_thm_verification_mesh();
realExpansion=build_thermo_interface_ring(S,3,p);
warm=solve_nonisothermal_ring(realExpansion,p,struct('cycles',2,'temperatureAmplitude_K',10, ...
    'saveFile',fullfile(folder,'case_mild_thermoelastic.mat')));
domainTestPassed=false;
try
    solve_nonisothermal_ring(realExpansion,p,struct('cycles',1,'temperatureAmplitude_K',20, ...
        'saveFile',fullfile(folder,'cold_thermoelastic_expected_domain_failure.mat')));
catch exception
    failurePath=fullfile(folder,'cold_thermoelastic_expected_domain_failure.mat.failure.mat');
    assert(isfile(failurePath),'ThermoInterface:FailureEvidence','Failure must preserve accepted history.');
    data=load(failurePath,'failure');
    domainTestPassed=strcmp(data.failure.domainDiagnostic.classification,'saturated_single_pressure_domain_exhausted');
    assert(domainTestPassed,'ThermoInterface:UnexpectedFailure','Unexpected failure: %s',exception.message);
end
assert(domainTestPassed,'ThermoInterface:DomainTest','Expected domain limit was not reached; review changed assumptions.');
limit=build_thermo_interface_ring(S,3,p,struct('thermalExpansion_Kinv',0));
reference=solve_nonisothermal_ring(limit,p,struct('temperatureAmplitude_K',20, ...
    'saveFile',fullfile(folder,'case_freezing_limit.mat')));
fine=solve_nonisothermal_ring(limit,p,struct('temperatureAmplitude_K',20,'dtDays',0.0625, ...
    'maximumPhaseIncrement',0.0125,'maximumPressureIncrement_Pa',5e3, ...
    'saveFile',fullfile(folder,'case_freezing_limit_dtHalf.mat')));
closed=solve_nonisothermal_ring(limit,p,struct('cycles',1,'temperatureAmplitude_K',20,'sealed',true, ...
    'saveFile',fullfile(folder,'case_sealed_limit.mat')));
report=struct('unitTestsPassed',sum([tests.Passed]),'domainLimitEvidencePreserved',domainTestPassed, ...
    'completedCases',{{warm.engine,reference.engine,fine.engine,closed.engine}}, ...
    'scope','Mixed-dimensional ring-interface verification only; surrounding soil and nonlinear contact are not integrated.', ...
    'notComplete20YearValidation',true);
save(fullfile(folder,'verification.mat'),'report','-v7');
fprintf('Completed verification cases; retain the cold thermoelastic domain failure separately.\n');
end
