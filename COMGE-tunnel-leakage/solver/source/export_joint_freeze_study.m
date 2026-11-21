function files = export_joint_freeze_study(study,folder)
%EXPORT_JOINT_FREEZE_STUDY Reproducible static figures and tables for v12.
%   STUDY.coldIce, .coldViscosity and .warmIce are complete twenty-year
%   RUN_TH_SEGMENT_YEAR20_JOINTFREEZE results. An optional .legacyConstant
%   supplies the previous constant-joint-conductance twenty-year result.
%   No field is invented if a legacy result lacks a diagnostic or runtime.
%
%   Leakage bins are the solver's conservative accepted-step integrals.
%   Daily joint diagnostics are time averages of accepted right-end states,
%   rather than reinterpolated leakage rates. Monthly leakage uses calendar
%   months; the twelve thermal snapshots use equal 365/12-day phases.
%   January 1 is phase day zero, the minimum of the supplied sinusoidal
%   tunnel-air boundary. Measured forcing, when supplied, takes precedence.
%
%   The radius panel is an illustrative constitutive-law comparison. It is
%   not a new coupled simulation or a measured gasket calibration. See
%   JOINT_FREEZE_CONDUCTANCE and its IAPWS/Gibbs-Thomson source references.
if nargin < 2 || isempty(folder)
    folder = fullfile(fileparts(mfilename('fullpath')),'joint_freeze_study_figures');
end
assert(isstruct(study) && isscalar(study), ...
    'export_joint_freeze_study:InvalidStudy','Study must be a scalar structure.');
required = {'coldIce','coldViscosity','warmIce'};
assert(all(isfield(study,required)), ...
    'export_joint_freeze_study:MissingMode', ...
    'Study needs coldIce, coldViscosity and warmIce results.');
modeNames = {'coldViscosity','coldIce','warmIce'};
modeLabels = {'Cold: viscosity','Cold: ice + viscosity','Warm: ice + viscosity'};
colors = [0.85 0.40 0.10;0.05 0.40 0.70;0.20 0.60 0.35];
if isfield(study,'legacyConstant') && ~isempty(study.legacyConstant)
    modeNames = [{'legacyConstant'},modeNames];
    modeLabels = [{'Cold: constant C (legacy)'},modeLabels];
    colors = [0.45 0.45 0.45;colors];
end
nm = numel(modeNames);
results = cell(nm,1);
daily = table(); monthly = table(); summary = table(); annual = table();
dailyByMode = cell(nm,1);
for mode = 1:nm
    r = study.(modeNames{mode});
    assert(isfield(r,'completedYears') && r.completedYears == 20 && ...
        isfield(r,'lastYear') && ~isempty(r.lastYear), ...
        'export_joint_freeze_study:IncompleteRun', ...
        'Mode %s has not completed its full twenty-year history.',modeNames{mode});
    results{mode} = r;
    out = r.lastYear; leak = out.leakage;
    assert(numel(leak.dailyTotalVolume_m3_per_m) == 365 && ...
        numel(leak.monthlyTotalVolume_m3_per_m) == 12, ...
        'export_joint_freeze_study:IncompleteBins', ...
        'Mode %s needs a full calendar year of daily and monthly bins.',modeNames{mode});
    dayEnd = leak.dailyBinEndDay(:) - 19*365;
    assert(max(abs(dayEnd-(1:365)')) < 1e-6, ...
        'export_joint_freeze_study:WrongYear','Expected year-20 daily bins.');
    day = dayEnd-0.5;
    q = 1000*leak.dailyTotalVolume_m3_per_m(:);
    [tmin,tmax,factor] = dailyFreezeDiagnostics(out,7);
    d = table(repmat(string(modeNames{mode}),365,1),day,q,tmin,tmax,factor, ...
        'VariableNames',{'mode','day_of_year_midpoint', ...
        'total_leakage_L_per_day_per_m','joint7_path_minT_dailyMean_C', ...
        'joint7_path_maxT_dailyMean_C','joint7_factor_dailyMean'});
    if mode == 1, daily = d; else, daily = [daily;d]; end %#ok<AGROW>
    dailyByMode{mode} = d;
    m = table(repmat(string(modeNames{mode}),12,1),(1:12)', ...
        1000*leak.monthlyTotalVolume_m3_per_m(:), ...
        'VariableNames',{'mode','calendar_month','total_leakage_L_per_m'});
    if mode == 1, monthly = m; else, monthly = [monthly;m]; end %#ok<AGROW>
    waterIndex = historyColumn(out,'waterError_m2');
    energyIndex = historyColumn(out,'energyError_J_per_m');
    runtime = NaN;
    if isfield(r,'elapsedWallTime_s'), runtime = r.elapsedWallTime_s; end
    rejected = NaN;
    if isfield(out,'solverWork') && isfield(out.solverWork,'rejectedSteps')
        rejected = out.solverWork.rejectedSteps;
    end
    a = annualVolumes(r,modeNames{mode});
    if mode == 1, annual = a; else, annual = [annual;a]; end %#ok<AGROW>
    acceptedTotal = NaN;
    if isfield(r,'columns') && any(strcmp(r.columns,'acceptedSteps'))
        index = find(strcmp(r.columns,'acceptedSteps'),1);
        acceptedTotal = sum(r.diagnostics(:,index));
    end
    s = table(string(modeNames{mode}),1000*leak.totalVolume_m3_per_m,runtime, ...
        max(abs(out.history(:,waterIndex))),max(abs(out.history(:,energyIndex))), ...
        size(out.history,1),rejected,acceptedTotal, ...
        'VariableNames',{'mode','year20_annual_leakage_L_per_m','runtime_s', ...
        'year20_max_abs_step_water_defect_m3_per_m', ...
        'year20_max_abs_step_energy_defect_J_per_m', ...
        'year20_accepted_steps','year20_rejected_steps','full20_accepted_steps'});
    if mode == 1, summary = s; else, summary = [summary;s]; end %#ok<AGROW>
end

cold = study.coldIce.lastYear;
nj = numel(cold.parameters.leakage.jointAnglesRad);
assert(nj >= 7 && isfield(cold.leakage,'jointThermalPaths') && ...
    isfield(cold.leakage,'jointFreeze'), ...
    'export_joint_freeze_study:MissingThermalPaths', ...
    'The coldIce result must contain joint-7 thermal paths and diagnostics.');
paths = cold.leakage.jointThermalPaths;
assert(numel(cold.snapshots) == 12, ...
    'export_joint_freeze_study:MissingSnapshots', ...
    'The mechanism figure requires twelve equal-phase snapshots.');
nq = paths.nq;
distance = ((1:nq)-0.5)/nq * ...
    (paths.radiusOuter_m(7)-paths.radiusInner_m(7));
profile = zeros(12,nq); snapshotDay = zeros(12,1);
for phase = 1:12
    allProfile = reshape(paths.interpolation*cold.snapshots(phase).T,nq,nj).';
    profile(phase,:) = allProfile(7,:);
    snapshotDay(phase) = cold.snapshots(phase).timeDay-19*365;
end
thermalPath = table(repelem((1:12)',nq),repelem(snapshotDay,nq), ...
    repmat((1:nq)',12,1),repmat(distance(:),12,1),reshape(profile.',[],1), ...
    'VariableNames',{'equal_phase_index','day_of_year','path_midpoint_index', ...
    'distance_from_outer_face_m','joint7_temperature_C'});

curveTemperature = linspace(-1,1,401)';
radiusValues = [0.1 1 10]*1e-6;
curveFactors = zeros(numel(curveTemperature),numel(radiusValues));
curveTable = table();
settings = cold.parameters.leakage.jointFreeze;
settings.enabled = true; settings.iceEnabled = true; settings.viscosityEnabled = true;
if ~isscalar(settings.logRadiusStd), settings.logRadiusStd = settings.logRadiusStd(7); end
for k = 1:numel(radiusValues)
    settings.flowWeightedRadius_m = radiusValues(k);
    curveFactors(:,k) = joint_freeze_conductance( ...
        repmat(curveTemperature,1,nq),settings);
    c = table(repmat(radiusValues(k),numel(curveTemperature),1), ...
        curveTemperature,curveFactors(:,k), ...
        'VariableNames',{'illustrative_flowWeightedRadius_m','uniform_path_T_C', ...
        'conductance_factor_relative_to_reference_T'});
    if k == 1, curveTable = c; else, curveTable = [curveTable;c]; end %#ok<AGROW>
end

if ~isfolder(folder), mkdir(folder); end
files = struct('dailyCSV',fullfile(folder,'year20_daily.csv'), ...
    'monthlyCSV',fullfile(folder,'year20_monthly.csv'), ...
    'summaryCSV',fullfile(folder,'year20_annual_summary.csv'), ...
    'annualHistoryCSV',fullfile(folder,'annual_leakage_year1_to20.csv'), ...
    'thermalPathCSV',fullfile(folder,'year20_joint7_thermal_path.csv'), ...
    'freezeLawCSV',fullfile(folder,'freeze_law_curves.csv'), ...
    'seasonalPNG',fullfile(folder,'figure1_seasonal_leakage.png'), ...
    'seasonalPDF',fullfile(folder,'figure1_seasonal_leakage.pdf'), ...
    'mechanismPNG',fullfile(folder,'figure2_freeze_thaw_mechanism.png'), ...
    'mechanismPDF',fullfile(folder,'figure2_freeze_thaw_mechanism.pdf'));
writetable(daily,files.dailyCSV,'Encoding','UTF-8');
writetable(monthly,files.monthlyCSV,'Encoding','UTF-8');
writetable(summary,files.summaryCSV,'Encoding','UTF-8');
writetable(annual,files.annualHistoryCSV,'Encoding','UTF-8');
writetable(thermalPath,files.thermalPathCSV,'Encoding','UTF-8');
writetable(curveTable,files.freezeLawCSV,'Encoding','UTF-8');

f1 = figure('Visible','off','Color','w','Position',[80 80 1350 880]);
clean1 = onCleanup(@() closeIfValid(f1));
layout1 = tiledlayout(f1,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout1); hold(ax,'on');
for mode = 1:nm
    d = dailyByMode{mode};
    plot(ax,d.day_of_year_midpoint,d.total_leakage_L_per_day_per_m, ...
        'Color',colors(mode,:),'LineWidth',1.6,'DisplayName',modeLabels{mode});
end
title(ax,'(a) Year-20 daily leakage');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'Leakage (L day^{-1} m^{-1})');
xlim(ax,[0 365]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);

coldIndex = find(strcmp(modeNames,'coldIce'),1);
dcold = dailyByMode{coldIndex};
ax = nexttile(layout1); hold(ax,'on');
plot(ax,dcold.day_of_year_midpoint,dcold.joint7_path_minT_dailyMean_C, ...
    'Color',[0.10 0.40 0.75],'LineWidth',1.6,'DisplayName','Joint-7 path minimum');
plot(ax,dcold.day_of_year_midpoint,dcold.joint7_path_maxT_dailyMean_C, ...
    'Color',[0.85 0.30 0.15],'LineWidth',1.6,'DisplayName','Joint-7 path maximum');
plot(ax,dcold.day_of_year_midpoint,tunnelAirAt(dcold.day_of_year_midpoint,cold.parameters), ...
    'k--','LineWidth',1.2,'DisplayName','Tunnel air');
yline(ax,cold.parameters.leakage.jointFreeze.bulkFreezingTemperature_C, ...
    ':','HandleVisibility','off');
title(ax,'(b) Cold case: joint-7 thermal state');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'Temperature (^{\circ}C)');
xlim(ax,[0 365]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);

ax = nexttile(layout1); hold(ax,'on');
plot(ax,dcold.day_of_year_midpoint,dcold.joint7_factor_dailyMean, ...
    'Color',colors(coldIndex,:),'LineWidth',1.8);
title(ax,'(c) Cold case: reversible joint conductance');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'C_{eff}/C_{reference} (daily mean)');
xlim(ax,[0 365]); ylim(ax,[0 max(1,1.05*max(dcold.joint7_factor_dailyMean))]); formatAxes(ax);

ax = nexttile(layout1); hold(ax,'on');
monthlyMatrix = zeros(12,nm);
for mode = 1:nm
    monthlyMatrix(:,mode) = 1000*results{mode}.lastYear.leakage.monthlyTotalVolume_m3_per_m(:);
end
bars = bar(ax,1:12,monthlyMatrix,'grouped');
for mode = 1:nm, bars(mode).FaceColor = colors(mode,:); end
title(ax,'(d) Year-20 calendar-month leakage');
xlabel(ax,'Calendar month'); ylabel(ax,'Monthly leakage (L m^{-1})');
xticks(ax,1:12); xticklabels(ax,monthLabels()); xlim(ax,[0.5 12.5]);
legend(ax,bars,modeLabels,'Location','best','FontSize',9); formatAxes(ax);
title(layout1,'Longitudinal-joint leakage under repeated freeze-thaw forcing: illustrative parameters', ...
    'FontName','Arial','FontSize',15,'FontWeight','normal');
exportgraphics(f1,files.seasonalPNG,'Resolution',300);
exportgraphics(f1,files.seasonalPDF,'ContentType','vector');
clear clean1

f2 = figure('Visible','off','Color','w','Position',[80 80 1350 880]);
clean2 = onCleanup(@() closeIfValid(f2));
layout2 = tiledlayout(f2,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout2);
imagesc(ax,distance,1:12,profile); set(ax,'YDir','normal');
bound = max(abs(profile(:))); if bound == 0, bound = 1; end
clim(ax,[-bound bound]); colormap(ax,temperatureMap(256));
cb = colorbar(ax); cb.Label.String = 'Temperature (^{\circ}C)';
title(ax,'(a) Joint-7 radial thermal path (cold case)');
xlabel(ax,'Distance from outer lining face (m; outside \rightarrow inside)');
xlim(ax,[0 paths.radiusOuter_m(7)-paths.radiusInner_m(7)]);
ylabel(ax,'Snapshot phase: 12 equal intervals'); yticks(ax,1:12);
formatAxes(ax); grid(ax,'off');

ax = nexttile(layout2); hold(ax,'on');
histDay = cold.history(:,historyColumn(cold,'day'))-19*365;
x = cold.leakage.jointFreeze.minTemperature_C(:,7);
y = cold.leakage.jointRate_m2_per_s(:,7)*86400*1000;
plot(ax,x,y,'Color',[0.75 0.75 0.75],'LineWidth',0.8,'HandleVisibility','off');
scatter(ax,x,y,24,histDay,'filled');
colormap(ax,parula(256)); clim(ax,[0 365]); cb = colorbar(ax);
cb.Label.String = 'Day of year (Jan 1 = 0)';
title(ax,'(b) Seasonal path: joint-7 leakage and coldest point');
xlabel(ax,'Path minimum temperature (^{\circ}C)'); ylabel(ax,'Joint leakage (L day^{-1} m^{-1})');
formatAxes(ax);

ax = nexttile(layout2); hold(ax,'on');
for mode = 1:nm
    a = annual(annual.mode == string(modeNames{mode}),:);
    plot(ax,a.year,a.annual_leakage_L_per_m,'-o','Color',colors(mode,:), ...
        'LineWidth',1.4,'MarkerSize',3,'DisplayName',modeLabels{mode});
end
title(ax,'(c) Annual leakage over the full twenty-year history');
xlabel(ax,'Simulation year'); ylabel(ax,'Annual leakage (L m^{-1})');
xlim(ax,[1 20]); xticks(ax,[1 5 10 15 20]);
legend(ax,'Location','best','FontSize',9); formatAxes(ax);

ax = nexttile(layout2); hold(ax,'on');
set(ax,'YScale','log');
radiusColors = [0.20 0.60 0.35;0.05 0.40 0.70;0.85 0.40 0.10];
for k = 1:numel(radiusValues)
    semilogy(ax,curveTemperature,curveFactors(:,k),'Color',radiusColors(k,:), ...
        'LineWidth',1.7,'DisplayName',[sprintf('r_Q = %g ',radiusValues(k)*1e6) '\mum']);
end
title(ax,{'(d) Illustrative microchannel constitutive curves', ...
    '(uniform temperature; not a coupled radius scan)'});
xlabel(ax,'Uniform path temperature (^{\circ}C)'); ylabel(ax,'C_{eff}/C_{reference}');
xlim(ax,[-1 1]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);
title(layout2,'Mechanism diagnostics: reversible blockage, thermal delay and transient history', ...
    'FontName','Arial','FontSize',15,'FontWeight','normal');
exportgraphics(f2,files.mechanismPNG,'Resolution',300);
exportgraphics(f2,files.mechanismPDF,'ContentType','vector');
clear clean2
end

function index = historyColumn(out,name)
assert(isfield(out,'historyColumns') && any(strcmp(out.historyColumns,name)), ...
    'export_joint_freeze_study:MissingHistoryColumn', ...
    'The result lacks history column %s.',name);
index = find(strcmp(out.historyColumns,name),1);
end

function [tmin,tmax,factor] = dailyFreezeDiagnostics(out,joint)
tmin = NaN(365,1); tmax = tmin; factor = tmin;
if ~isfield(out.leakage,'jointFreeze'), return; end
jf = out.leakage.jointFreeze;
if size(jf.factor,2) < joint, return; end
time = out.history(:,historyColumn(out,'day'));
dt = out.history(:,historyColumn(out,'dt_day'));
assert(size(jf.factor,1) == numel(time), ...
    'export_joint_freeze_study:InconsistentDiagnostics', ...
    'Joint-freezing history must match accepted steps.');
values = [jf.minTemperature_C(:,joint),jf.maxTemperature_C(:,joint),jf.factor(:,joint)];
edges = 19*365+(0:365)';
assert(max(abs(diff([edges(1);time])-dt)) < 1e-6, ...
    'export_joint_freeze_study:HistoryGap','Annual diagnostic history has a time gap.');
cumulative = [zeros(1,3);cumsum(values.*dt,1)];
dailyMeans = diff(interp1([edges(1);time],cumulative,edges,'linear'),1,1);
tmin = dailyMeans(:,1); tmax = dailyMeans(:,2); factor = dailyMeans(:,3);
end

function a = annualVolumes(result,mode)
assert(isfield(result,'diagnostics') && isfield(result,'columns') && ...
    all(ismember({'year','leak_m3_per_m'},result.columns)), ...
    'export_joint_freeze_study:MissingAnnualHistory', ...
    'Mode %s lacks the full annual diagnostic history.',mode);
yearIndex = find(strcmp(result.columns,'year'),1);
volumeIndex = find(strcmp(result.columns,'leak_m3_per_m'),1);
year = result.diagnostics(:,yearIndex);
assert(isequal(year(:),(1:20)'), ...
    'export_joint_freeze_study:IncompleteAnnualHistory', ...
    'Mode %s must retain annual volumes for years 1 through 20.',mode);
a = table(repmat(string(mode),20,1),year,1000*result.diagnostics(:,volumeIndex), ...
    'VariableNames',{'mode','year','annual_leakage_L_per_m'});
end

function temperature = tunnelAirAt(day,p)
phase = mod(day,p.yearDays);
temperature = (p.tunnelMax+p.tunnelMin)/2 - ...
    (p.tunnelMax-p.tunnelMin)/2*cos(2*pi*phase/p.yearDays);
if isfield(p,'forcing') && isfield(p.forcing,'tunnelAirC') && ~isempty(p.forcing.tunnelAirC)
    temperature = interp1(p.forcing.day,p.forcing.tunnelAirC,phase,'linear');
end
end

function formatAxes(ax)
set(ax,'FontName','Arial','FontSize',11,'LineWidth',0.8,'Box','on', ...
    'TickDir','out','XMinorTick','off','YMinorTick','off');
grid(ax,'on'); ax.GridAlpha = 0.15;
end

function labels = monthLabels()
labels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
end

function colors = temperatureMap(n)
% Symmetric blue-white-red temperature map; zero always uses its midpoint.
anchors = [0.10 0.30 0.70;0.45 0.70 0.90;1 1 1;0.95 0.65 0.45;0.70 0.10 0.10];
colors = interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n),'linear');
end

function closeIfValid(fig)
if isgraphics(fig), close(fig); end
end
