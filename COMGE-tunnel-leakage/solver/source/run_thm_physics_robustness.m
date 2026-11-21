function report = run_thm_physics_robustness(user)
%RUN_THM_PHYSICS_ROBUSTNESS Independent histories; illustrative uncertainty.
% Default two-year qualification of actual R7 THM, local freeze/thaw memory,
% conductance/radius geometry and 16/32-point joint path integration.
% This is not field calibration. Parameter ranges are exploratory scenarios.
% Full histories can be requested with maxYears=20 in a fresh directory.
arguments
    user (1,1) struct = struct()
end
root = fileparts(mfilename('fullpath'));
config = struct('maxYears',2,'resume',false, ...
    'outputFolder',fullfile(root,'physical_robustness_20261003','R7_two_years'));
keys = fieldnames(user);
for index = 1:numel(keys)
    assert(isfield(config,keys{index}), 'THMRobustness:Option','Unknown option: %s.',keys{index});
    config.(keys{index}) = user.(keys{index});
end
validateattributes(config.maxYears,{'numeric'},{'scalar','integer','>=',2,'<=',20});
folder = char(config.outputFolder);
if ~java.io.File(folder).isAbsolute()
    folder = fullfile(pwd,folder);
end
folder = char(java.io.File(folder).getCanonicalPath());
if ~isfolder(folder)
    mkdir(folder);
end
config.outputFolder = folder;
file = fullfile(folder,'robustness_configuration.mat');
if isfile(file)
    old = load(file,'configuration');
    assert(config.resume && isequaln(rmfield(old.configuration,'resume'),rmfield(config,'resume')), ...
        'THMRobustness:Resume','Existing campaign requires identical configuration and resume=true.');
else
    configuration = config;
    save(file,'configuration','-v7');
end
% Paired nominal/control first, then integration refinement and uncertainty.
name = ["control","nominal","path32","rate_low","rate_high", ...
    "amplitude_low","amplitude_high","radius_small","radius_large"];
rate = [0 .05 .05 .02 .1 .05 .05 .05 .05];
alpha = [.5 .5 .5 .5 .5 .25 1 .5 .5];
radius = [1 1 1 1 1 1 1 .3 3]*1e-6;
pathCount = [16 16 32 16 16 16 16 16 16];
numberOfCases = numel(name);
results = cell(numberOfCases,1);
rows = zeros(numberOfCases,10);
for index = 1:numberOfCases
    caseFolder = fullfile(folder,char(name(index)));
    fprintf('ROBUSTNESS case %d/%d: %s, independent %d-year history.\n',index,numberOfCases,name(index),config.maxYears);
    progress = struct('caseIndex',index,'caseName',name(index),'completedCases',index-1, ...
        'status','running','config',config);
    save(fullfile(folder,'robustness_progress.mat'),'progress','-v7');
    conductance = zeros(1,8);
    conductance(7) = 2e-10*(radius(index)/1e-6)^2;
    % Radius scenarios keep the connected volume unchanged, rather than
    % independently changing C and radius and inadvertently altering storage.
    options = struct('parameterPack','kanagawa_th','maxYears',config.maxYears, ...
        'doseRates',rate(index),'cycleSettings',struct('radiusGrowthAtUnitDamage',alpha(index)), ...
        'model',struct('leakage',struct('jointConductance',conductance,'jointFreeze',struct( ...
        'quadraturePoints',pathCount(index),'flowWeightedRadius_m',radius(index)))), ...
        'outputFolder',caseFolder,'useFastNumerics',true, ...
        'resume',config.resume && isfile(fullfile(caseFolder,'cycle_study_configuration.mat')));
    run_thm_literature_parameter_study(options);
    raw = load(fullfile(caseFolder,'rate_01',sprintf('year_%02d_annual_results.mat',config.maxYears)),'out');
    out = raw.out;
    memory = out.jointCycleDamage.finalMemory;
    geometry = out.leakage.channelGeometry;
    rows(index,:) = [rate(index),alpha(index),radius(index)*1e6,pathCount(index), ...
        1000*out.leakage.totalVolume_m3_per_m,out.jointCycleDamage.appliedSeverity(7), ...
        mean(memory.cycleCount(7,:)),geometry.equivalentVoidFraction(7), ...
        out.solverWork.acceptedSteps,out.solverWork.rejectedSteps];
    results{index} = out;
    summary = array2table(rows(1:index,:), 'VariableNames', ...
        {'doseRate','alpha','radius_um','pathPoints','annualLeak_L_m','J7Severity', ...
        'J7MeanCycles','equivalentVoidFraction','acceptedSteps','rejectedSteps'});
    summary.caseName = name(1:index)';
    writetable(summary,fullfile(folder,'robustness_summary.csv'));
end
coarse = results{2};
fine = results{3};
pathAudit = comparePaths(coarse,fine);
summary.relativeToControl_percent = 100*(summary.annualLeak_L_m/summary.annualLeak_L_m(1)-1);
report = struct('passed',pathAudit.passed,'configuration',config,'summary',summary, ...
    'pathAudit',pathAudit,'parameterRangesAreCalibrated',false, ...
    'isFull20HistoryStudy',config.maxYears==20, ...
    'scope','Independent numerical scenario histories with literature soil; physical admission and path refinement; no empirical joint damage or field validation.');
save(fullfile(folder,'robustness_report.mat'),'report','-v7');
writetable(summary,fullfile(folder,'robustness_summary.csv'));
progress.status = 'completed';
progress.completedCases = numberOfCases;
progress.pathGatePassed = pathAudit.passed;
save(fullfile(folder,'robustness_progress.mat'),'progress','-v7');
disp(summary);
assert(report.passed,'THMRobustness:PathRefinement', ...
    'Joint path refinement gate failed. Evidence is retained; investigate integration before further prediction.');
end

function report = comparePaths(coarse,fine)
coarseParameters = coarse.parameters;
fineParameters = fine.parameters;
coarseParameters.leakage.jointFreeze.quadraturePoints = fineParameters.leakage.jointFreeze.quadraturePoints;
assert(isequaln(thm_constitutive_signature(coarseParameters), ...
    thm_constitutive_signature(fineParameters)) ...
    &&strcmp(coarseParameters.meshFile,fineParameters.meshFile) ...
    &&isequal(coarseParameters.integration,fineParameters.integration) ...
    &&coarse.timeDay==fine.timeDay, 'THMRobustness:PathPhysics', ...
    'Path refinement must compare the same physical case, mesh and complete duration.');
annualError = abs(coarse.leakage.totalVolume_m3_per_m/fine.leakage.totalVolume_m3_per_m-1);
monthError = max(abs(coarse.leakage.monthlyTotalVolume_m3_per_m ...
    -fine.leakage.monthlyTotalVolume_m3_per_m))/fine.leakage.totalVolume_m3_per_m;
meanCountDifference = max(abs(mean(coarse.jointCycleDamage.finalMemory.cycleCount,2) ...
    -mean(fine.jointCycleDamage.finalMemory.cycleCount,2)));
damageDifference = max(abs(coarse.jointCycleDamage.appliedSeverity-fine.jointCycleDamage.appliedSeverity));
report = struct('passed',annualError<=.01 && monthError<=.01 && meanCountDifference<1e-12, ...
    'annualRelativeFineError',annualError,'maximumMonthlyDifferenceOverFineAnnual',monthError, ...
    'maximumMeanCycleCountDifference',meanCountDifference,'maximumSeverityDifference',damageDifference, ...
    'scope','Complete same-duration 16/32 path-point histories; fixed mesh, physical parameters and tolerances. Not a spatial mesh convergence check.');
end
