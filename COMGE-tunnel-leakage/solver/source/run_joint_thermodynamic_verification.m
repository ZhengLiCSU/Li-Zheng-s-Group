function report = run_joint_thermodynamic_verification(folder)
%RUN_JOINT_THERMODYNAMIC_VERIFICATION Independent branch, never old THM evidence.
arguments
    folder (1,1) string = fullfile(fileparts(mfilename('fullpath')), ...
        'thermodynamic_interface_20261003_final_v2')
end
assert(~isfile(fullfile(folder,'verification.mat')),'JointThermo:ExistingVerification', ...
    'Preserve the existing verification; choose a new output folder.');
if ~isfolder(folder), mkdir(folder); end
sourceFolder = fullfile(folder,'source_at_start');
if ~isfolder(sourceFolder), mkdir(sourceFolder); end
sources = {'joint_thermodynamic_parameters.m','joint_helmholtz_state.m', ...
    'joint_interface_residual.m','joint_interface_primary_residual.m', ...
    'joint_thermodynamic_step.m','joint_closed_equilibrium.m', ...
    'joint_adiabatic_equilibrium.m','solve_thermodynamic_joint.m', ...
    'run_joint_thermodynamic_verification.m','TestJointThermodynamics.m'};
root = fileparts(mfilename('fullpath'));
for k = 1:numel(sources), copyfile(fullfile(root,sources{k}),fullfile(sourceFolder,sources{k})); end
testResults = runtests(fullfile(root,'TestJointThermodynamics.m'));
testNames = string({testResults.Name}); testPassed = [testResults.Passed];
assert(all(testPassed),'JointThermo:FailedTests','Do not run the study with failed kernel tests.');
names = ["reference","dtHalf","dtQuarter","space16","space32","space64","widerSlit","softWalls","noFreezing"];
cases = cell(numel(names),1);
base = struct('cells',16,'duration_s',2*365*86400,'dt_s',86400, ...
    'initialTemperature_K',285.15,'verbose',false);
commonTimeGrid = [];
for k = 1:numel(names)
    options = base;
    switch names(k)
        case "dtHalf", options.dt_s = 43200;
        case "dtQuarter", options.dt_s = 21600;
        case "space16", options.timeGrid_s = commonTimeGrid; options.dt_s = 21600;
        case "space32", options.cells = 32; options.timeGrid_s = commonTimeGrid; options.dt_s = 21600;
        case "space64", options.cells = 64; options.timeGrid_s = commonTimeGrid; options.dt_s = 21600;
        case "widerSlit", options.parameters = struct('aperture0_m',6e-6);
        case "softWalls", options.parameters = struct('normalStiffness_Pam',1e11);
        case "noFreezing", options.temperatureFunction = @(t)277.15;
    end
    fprintf('Independent thermostated slit case: %s\n',names(k));
    out = solve_thermodynamic_joint(options);
    if names(k)=="dtQuarter", commonTimeGrid = out.time_s; end
    save(fullfile(folder,names(k)+".mat"),'out','-v7');
    cases{k} = summarize(out);
end
fprintf('Independent ideal slit: 20 full temperature cycles (NOT full tunnel THM).\n');
options = base; options.duration_s = 20*365*86400; options.dt_s = 43200;
out = solve_thermodynamic_joint(options);
save(fullfile(folder,'twenty_cycles.mat'),'out','-v7');
twenty = summarize(out);
p = joint_thermodynamic_parameters();
[closedCold,closedReport] = joint_closed_equilibrium(0.005,272.15,p);
[adiabatic,adiabaticReport] = joint_adiabatic_equilibrium([0,0.005,0,272.15],p);
rigid = joint_thermodynamic_parameters(struct('normalStiffness_Pam',1e20));
[rigidCold,rigidReport] = joint_closed_equilibrium(0.005,272.15,rigid);
reference = cases{1}; fine = cases{3};
timeError = compare(cases{1},cases{2}); finerError = compare(cases{2},fine);
spaceError = compare(cases{4},cases{5}); finerSpaceError = compare(cases{5},cases{6});
report = struct('scope','Thermostated ideal slit; not full tunnel THM, not experimental gasket validation', ...
    'testNames',testNames,'testPassed',testPassed,'names',names,'cases',{cases}, ...
    'twentyCycles',twenty,'timeError',timeError,'finerTimeError',finerError, ...
    'spaceError',spaceError,'finerSpaceError',finerSpaceError,'commonSpatialTime_s',commonTimeGrid, ...
    'closedCold',closedCold,'closedPressure_Pa',closedReport.pressure_Pa, ...
    'adiabatic',adiabatic,'adiabaticAudit',adiabaticReport, ...
    'rigidCold',rigidCold,'rigidPressure_Pa',rigidReport.pressure_Pa, ...
    'referenceTargetYearQ_Lm',reference.annualVolume_Lm(2), ...
    'fineTargetYearQ_Lm',fine.annualVolume_Lm(2));
report.passed = all(testPassed) && timeError.annualRelative<0.01 && ...
    timeError.maxBinRelativeToAnnual<0.01 && finerError.annualRelative<0.01 && ...
    spaceError.annualRelative<0.01 && finerSpaceError.annualRelative<0.01;
save(fullfile(folder,'verification.mat'),'report','-v7');
fprintf('Slit verification passed=%d; Y2 %.12g L/m; dt error %.6g; mesh error %.6g.\n', ...
    report.passed,reference.annualVolume_Lm(2),timeError.annualRelative,spaceError.annualRelative);
end

function result = summarize(out)
year_s = 365*86400; endTime = out.time_s(2:end); startTime = out.time_s(1:end-1);
years = round(out.time_s(end)/year_s);
annual = zeros(years,1); bins = zeros(years,12);
% Exact overlap integration, independent of solver calendar alignment.
for year = 1:years
    for bin = 1:12
        low = ((year-1)+(bin-1)/12)*year_s;
        high = ((year-1)+bin/12)*year_s;
        overlap = max(0,min(endTime,high)-max(startTime,low));
        bins(year,bin) = 1000*sum(overlap.*out.outletReferenceFlux_m2s);
    end
    annual(year) = sum(bins(year,:));
end
result = struct('annualVolume_Lm',annual,'equalCalendarBins_Lm',bins, ...
    'massLedger_kgm',out.massLedger_kgm, ...
    'maxLocalMassDefect_kgm',max(out.audits(:,3)), ...
    'maxEnergyIdentityDefect_Jm',max(abs(out.audits(:,4))), ...
    'minimumDarcyDissipation_Jm',min(out.audits(:,5)), ...
    'minimumPhaseDissipation_Jm',min(out.audits(:,6)), ...
    'minimumBackwardEulerRemainder_Jm',min(out.audits(:,7)), ...
    'minimumEntropyProduction_JmK',min(out.audits(:,8)), ...
    'maxPressure_Pa',max(out.pressureHistory_Pa,[],'all'),'elapsed_s',out.elapsed_s, ...
    'steps',numel(out.dt_s),'rejectedSteps',out.rejectedSteps,'dtMin_s',min(out.dt_s), ...
    'cells',out.options.cells,'dtMax_s',out.options.dt_s);
end

function error = compare(coarse,fine)
error = struct('annualRelative',abs(coarse.annualVolume_Lm(2)-fine.annualVolume_Lm(2))/ ...
    fine.annualVolume_Lm(2),'maxBinRelativeToAnnual', ...
    max(abs(coarse.equalCalendarBins_Lm(2,:)-fine.equalCalendarBins_Lm(2,:)))/fine.annualVolume_Lm(2));
end
