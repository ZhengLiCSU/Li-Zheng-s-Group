function files=export_joint_supply_study(study,folder)
%EXPORT_JOINT_SUPPLY_STUDY Audited full-history leakage approximations.
% files=export_joint_supply_study loads all five LAUNCH_JOINT_SUPPLY_STUDY
% results. STUDY may instead be a study-folder path or a scalar struct with
% fields base, transition, supply, supplySeasonal, independentNoDrain.
% Only completed v14 twenty-year runs are accepted.
%
% The same-temperature approximation uses FULL TH joint temperatures, so
% its error isolates omitted supply head loss but does not remove thermal
% simulation. The independent thermal approximation uses the matched
% no-extraction TH history only for supplySeasonal, never for another soil
% permeability or boundary scenario. Neither approximation is a separate
% quasi-steady hydraulic solution, nor an estimate of hydraulic storage's
% isolated effect. Figures show reversible annual forcing and transient
% histories, not accumulating freeze-thaw damage or calibrated site data.
root=fileparts(mfilename('fullpath')); addpath(root);
defaultStudy=fullfile(root,'joint_supply_study_20260930');
if nargin<1 || isempty(study), study=defaultStudy; end
if ischar(study)||isstring(study)
 studyFolder=char(study); study=loadStudy(studyFolder);
else
 studyFolder=defaultStudy;
end
if nargin<2 || isempty(folder), folder=fullfile(studyFolder,'exports'); end
names={'base','transition','supply','supplySeasonal','independentNoDrain'};
assert(isstruct(study)&&isscalar(study)&&all(isfield(study,names)), ...
 'All five paired study cases are required.');
if ~isfolder(folder), mkdir(folder); end
series=cell(numel(names),1); summary=table(); daily=table(); monthly=table();
annual=table(); accepted=table();
background=study.independentNoDrain;
assert(~background.lastYear.parameters.tunnelDrain, ...
 'Independent thermal input must come from a no-extraction TH model.');
for k=1:numel(names)
 r=study.(names{k});
 assert(isfield(r,'completedYears')&&r.completedYears==20&& ...
  strcmp(r.lastYear.version,'TH_v14_20260930_joint_supply'), ...
  'Case %s requires a completed v14 twenty-year history.',names{k});
 s=buildSeries(r,names{k}); series{k}=s;
 daily=appendTable(daily,s.daily);
 monthly=appendTable(monthly,s.monthly);
 accepted=appendTable(accepted,s.accepted);
 h=r.lastYear.history; out=r.lastYear;
 water=h(:,historyColumn(out,'waterError_m2'));
 energy=h(:,historyColumn(out,'energyError_J_per_m'));
 q=s.Qactual; qSameStatic=s.QsameStatic; qSameSeasonal=s.QsameSeasonal;
 qIndependent=NaN; errIndependent=NaN; l1Independent=NaN;
 if strcmp(names{k},'supplySeasonal')
  validateBackground(r,background);
  bg=buildSeries(background,'independentNoDrain');
  qIndependent=bg.QsameSeasonal;
  errIndependent=relative(qIndependent-q,qIndependent);
  l1Independent=relative(sum(abs(bg.daily.sameT_seasonal_L_per_day_per_m- ...
   s.daily.actual_L_per_day_per_m)),q);
 end
 errStatic=NaN; errSeasonal=NaN; l1Static=NaN; l1Seasonal=NaN; ratio=NaN;
 if out.parameters.tunnelDrain
  errStatic=relative(qSameStatic-q,qSameStatic);
  errSeasonal=relative(qSameSeasonal-q,qSameSeasonal);
  ratio=relative(q,qSameSeasonal);
  l1Static=relative(sum(abs(s.daily.sameT_static_L_per_day_per_m- ...
   s.daily.actual_L_per_day_per_m)),q);
  l1Seasonal=relative(sum(abs(s.daily.sameT_seasonal_L_per_day_per_m- ...
   s.daily.actual_L_per_day_per_m)),q);
 end
 open=s.factor>1e-8;
 heffMin=NaN; heffMax=NaN;
 if any(open)&&out.parameters.tunnelDrain
  heffMin=min(s.actualHead(open)); heffMax=max(s.actualHead(open));
 end
 [peak,peakDay]=max(s.daily.actual_L_per_day_per_m);
 if peak<=1e-12, peakDay=NaN; end
 isFrozenDay=s.daily.factor_mean<=1e-8;
 nj=numel(out.parameters.leakage.jointAnglesRad);
 C=expandJoint(out.parameters.leakage.jointConductance,nj);
 runtime=NaN; if isfield(r,'elapsedWallTime_s'), runtime=r.elapsedWallTime_s; end
 periodic=false; if isfield(r,'periodicConverged'), periodic=r.periodicConverged; end
 row=table(string(names{k}),out.parameters.Ksat,C(7),out.parameters.tunnelDrain, ...
  q,s.QconstantStatic,qSameStatic,qSameSeasonal,qIndependent, ...
  100*errStatic,100*errSeasonal,100*errIndependent, ...
  100*l1Static,100*l1Seasonal,100*l1Independent,ratio, ...
  sum(s.Ceff.*s.dt_s),heffMin,heffMax,sum(isFrozenDay)/365,peak,peakDay, ...
  runtime,r.completedYears,periodic,size(out.pAll,1),numel(out.soilNodes), ...
  size(out.mesh.tSoil,1),size(h,1),out.solverWork.rejectedSteps, ...
  min(h(:,historyColumn(out,'dt_day'))),max(h(:,historyColumn(out,'scaledResidual'))), ...
  sum(water),sum(abs(water)),max(abs(water)),sum(energy),sum(abs(energy)), ...
  max(abs(energy)),max(out.leakage.hydraulics.negativeTransEdges), ...
  'VariableNames',{'caseName','Ksat_m_per_s','joint7_Cref_m_per_s','tunnelExtraction', ...
  'Q20_actual_L_per_m','Q20_constantC_static_L_per_m','Q20_sameT_static_L_per_m', ...
  'Q20_sameT_seasonal_L_per_m','Q20_independentT_seasonal_L_per_m', ...
  'static_overprediction_percent','seasonal_overprediction_percent', ...
  'independent_overprediction_percent','dailyL1_static_percentOfActual', ...
  'dailyL1_seasonal_percentOfActual','dailyL1_independent_percentOfActual', ...
  'actual_over_sameT_seasonal','Ceff_timeIntegral_m','openPeriod_head_min_m', ...
  'openPeriod_head_max_m','nearBlockedDailyFraction_factorBelow1em8', ...
  'peakDailyActual_L_per_day_per_m','peakDailyActual_binEndDay', ...
  'lastInvocation_runtime_s','completedYears','periodicConverged', ...
  'thermalNodes','soilHeadNodes','soilTriangles','year20_acceptedSteps', ...
  'year20_rejectedSteps','year20_minStep_day','year20_maxScaledResidual', ...
  'year20_signedWaterDefect_m3_per_m','year20_sumAbsWaterDefect_m3_per_m', ...
  'year20_maxAbsStepWaterDefect_m3_per_m','year20_signedEnergyDefect_J_per_m', ...
  'year20_sumAbsEnergyDefect_J_per_m','year20_maxAbsStepEnergyDefect_J_per_m', ...
  'year20_maxNegativeTransEdges'});
 summary=appendTable(summary,row);
 yy=r.diagnostics(:,find(strcmp(r.columns,'year'),1));
 qq=1000*r.diagnostics(:,find(strcmp(r.columns,'leak_m3_per_m'),1));
 assert(isequal(yy(:),(1:20)'), ...
  'Case %s must retain the entire annual history.',names{k});
 annual=appendTable(annual,table(repmat(string(names{k}),20,1),yy,qq, ...
  'VariableNames',{'caseName','year','annualActual_L_per_m'}));
end
files=struct('summaryCSV',fullfile(folder,'year20_supply_summary.csv'), ...
 'dailyCSV',fullfile(folder,'year20_supply_daily.csv'), ...
 'monthlyCSV',fullfile(folder,'year20_supply_monthly.csv'), ...
 'annualHistoryCSV',fullfile(folder,'annual_actual_year1_to20.csv'), ...
 'acceptedStepCSV',fullfile(folder,'year20_accepted_supply_diagnostics.csv'), ...
 'overviewPNG',fullfile(folder,'figure1_supply_regimes.png'), ...
 'overviewPDF',fullfile(folder,'figure1_supply_regimes.pdf'), ...
 'seasonalPNG',fullfile(folder,'figure2_seasonal_supply.png'), ...
 'seasonalPDF',fullfile(folder,'figure2_seasonal_supply.pdf'), ...
 'notesTXT',fullfile(folder,'READ_FIGURES.txt'));
writetable(summary,files.summaryCSV,'Encoding','UTF-8');
writetable(daily,files.dailyCSV,'Encoding','UTF-8');
writetable(monthly,files.monthlyCSV,'Encoding','UTF-8');
writetable(annual,files.annualHistoryCSV,'Encoding','UTF-8');
writetable(accepted,files.acceptedStepCSV,'Encoding','UTF-8');
makeOverview(series,summary,annual,files);
makeSeasonal(series{4},series{5},files);
fid=fopen(files.notesTXT,'w'); assert(fid>=0,'Cannot write figure notes.');
cleanFile=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,['All cases have independent full 20-year histories on a level-1 refined active-joint mesh.\n' ...
 'Parameters and groundwater oscillation are uncalibrated scenarios, not a site prediction.\n' ...
 'sameT means temperature diagnostics from the actual coupled TH solution.\n' ...
 'independentT uses the no-extraction background TH case with the SAME soil and groundwater forcing.\n' ...
 'The background case is matched only to supplySeasonal, not to the other cases.\n' ...
 'Neither approximation solves the separate quasi-steady hydraulic problem.\n' ...
 'Annual leakage evolution represents transient fields and reversible icing, not cumulative damage.\n' ...
 'Daily and monthly rates are conservative accepted-step bin integrals. Diagnostic means are time weighted.\n' ...
 'Near-blocked days use factor <= 1e-8; this is a reporting threshold, not a measured gasket closure.\n' ...
 'lastInvocation_runtime_s is full-history wall time only when the run did not resume.\n' ...
 'Three constant-water-table points form a one-dimensional permeability scan; they are not a two-dimensional phase map.\n']);
end

function study=loadStudy(folder)
names={'base','transition','supply','supplySeasonal','independentNoDrain'};
study=struct();
for k=1:numel(names)
 file=fullfile(folder,names{k},'TH_segment_spinup_results.mat');
 assert(isfile(file),'Missing completed case file: %s.',file);
 z=load(file,'result'); study.(names{k})=z.result;
end
end

function s=buildSeries(result,name)
out=result.lastYear; p=out.parameters; L=p.leakage;
nj=numel(L.jointAnglesRad);
assert(nj==8&&all(expandJoint(L.jointEntryHead,nj)==0), ...
 'This paired study requires the eight-joint, zero-entry-head scenario.');
C=expandJoint(L.jointConductance,nj);
assert(nnz(C)>0&&all(find(C>0)==7), ...
 'This paired study requires joint 7 as its single conducting joint.');
time=out.history(:,historyColumn(out,'day'));
dt=out.history(:,historyColumn(out,'dt_day'))*86400;
start=19*365;
assert(max(abs(diff([start;time])*86400-dt))<1e-5&& ...
 abs(time(end)-20*365)<1e-7,'Year-20 accepted time history must be complete.');
f=out.leakage.jointFreeze.factor(:,7);
hyd=out.leakage.hydraulics;
assert(size(hyd.history,1)==numel(time)&&size(hyd.drivingHead_m,2)==nj&& ...
 all(isfinite(hyd.history(:))),'Hydraulic diagnostic history is incomplete.');
Ceff=C(7)*f;
q=out.leakage.jointRate_m2_per_s(:,7);
qRef=hyd.unrestrictedRate_m2_per_s(:,7);
staticHead=hyd.referenceHead_m(:,7)+hyd.waterTableDepth_m-p.waterTableDepth;
qStatic=Ceff.*max(staticHead,0);
qConstant=C(7)*max(staticHead,0);
qCheck=Ceff.*max(hyd.referenceHead_m(:,7),0);
assert(max(abs(qCheck-qRef))<1e-17+1e-8*max(abs(qRef)), ...
 'Unrestricted leakage diagnostics and the defined reference head disagree.');
actualHead=hyd.drivingHead_m(:,7);
phase=time-start;
values=[q,qStatic,qRef,qConstant];
dailyVolume=binAccepted(time,dt,values,start+(0:365)')*1000;
monthEdges=start+[0 31 59 90 120 151 181 212 243 273 304 334 365]';
monthlyVolume=binAccepted(time,dt,values,monthEdges)*1000;
assert(max(abs(dailyVolume(:,1)-1000*out.leakage.dailyTotalVolume_m3_per_m(:)))<1e-7, ...
 'Daily bins must reproduce the conservative solver leakage.');
assert(abs(sum(dailyVolume(:,1))-1000*out.leakage.totalVolume_m3_per_m)<1e-7, ...
 'Annual leakage does not equal the daily integral.');
diagnostics=[f,Ceff,actualHead,hyd.referenceHead_m(:,7),staticHead, ...
 out.leakage.jointFreeze.minTemperature_C(:,7), ...
 out.leakage.jointFreeze.maxTemperature_C(:,7),hyd.waterTableDepth_m];
dailyMeans=binAccepted(time,dt,diagnostics,start+(0:365)')/86400;
daily=table(repmat(string(name),365,1),(1:365)', ...
 dailyVolume(:,1),dailyVolume(:,2),dailyVolume(:,3),dailyVolume(:,4), ...
 dailyMeans(:,1),dailyMeans(:,2),dailyMeans(:,3),dailyMeans(:,4), ...
 dailyMeans(:,5),dailyMeans(:,6),dailyMeans(:,7),dailyMeans(:,8), ...
 'VariableNames',{'caseName','day_of_year_binEnd','actual_L_per_day_per_m', ...
 'sameT_static_L_per_day_per_m','sameT_seasonal_L_per_day_per_m', ...
 'constantC_static_L_per_day_per_m','factor_mean','Ceff_mean_m_per_s', ...
 'actualDrivingHead_mean_m','seasonalReferenceHead_mean_m','staticReferenceHead_mean_m', ...
 'pathMinT_mean_C','pathMaxT_mean_C','waterTableDepth_mean_m'});
monthly=table(repmat(string(name),12,1),(1:12)',monthlyVolume(:,1), ...
 monthlyVolume(:,2),monthlyVolume(:,3),monthlyVolume(:,4), ...
 'VariableNames',{'caseName','month','actual_L_per_m','sameT_static_L_per_m', ...
 'sameT_seasonal_L_per_m','constantC_static_L_per_m'});
accepted=table(repmat(string(name),numel(time),1),phase,dt/86400,q,qStatic,qRef, ...
 f,Ceff,actualHead,hyd.referenceHead_m(:,7),staticHead,hyd.waterTableDepth_m, ...
 hyd.minTrans_m_per_s,hyd.negativeTransEdges, ...
 'VariableNames',{'caseName','day_of_year_end','dt_day','actual_m2_per_s', ...
 'sameT_static_m2_per_s','sameT_seasonal_m2_per_s','factor','Ceff_m_per_s', ...
 'actualDrivingHead_m','seasonalReferenceHead_m','staticReferenceHead_m', ...
 'waterTableDepth_m','minTrans_m_per_s','negativeTransEdges'});
s=struct('name',name,'daily',daily,'monthly',monthly,'accepted',accepted, ...
 'dt_s',dt,'factor',f,'Ceff',Ceff,'actualHead',actualHead, ...
 'Qactual',1000*sum(q.*dt),'QsameStatic',1000*sum(qStatic.*dt), ...
 'QsameSeasonal',1000*sum(qRef.*dt),'QconstantStatic',1000*sum(qConstant.*dt));
end

function validateBackground(actual,background)
a=actual.lastYear.parameters; b=background.lastYear.parameters;
assert(a.Ksat==b.Ksat&&a.waterTableDepth==b.waterTableDepth&& ...
 isequaln(a.forcing,b.forcing)&&isequaln(a.soilViscosity,b.soilViscosity)&& ...
 a.topMin==b.topMin&&a.topMax==b.topMax&&a.tunnelMin==b.tunnelMin&& ...
 a.tunnelMax==b.tunnelMax&&isequaln(a.leakage,b.leakage)&& ...
 isequaln(actual.lastYear.pAll,background.lastYear.pAll), ...
 'The independent thermal approximation must use matched materials, forcing and mesh.');
assert(a.tunnelDrain&&~b.tunnelDrain,'Extraction flags are not a matched pair.');
end

function volume=binAccepted(time,dt_s,rates,edges)
% Each accepted rate is held constant over its own implicit Euler interval.
first=time(1)-dt_s(1)/86400;
cumulative=[zeros(1,size(rates,2));cumsum(rates.*dt_s,1)];
volume=diff(interp1([first;time],cumulative,edges,'linear'),1,1);
assert(all(isfinite(volume(:))),'Conservative binning encountered an incomplete interval.');
end

function index=historyColumn(out,name)
index=find(strcmp(out.historyColumns,name),1);
assert(~isempty(index),'Missing history column %s.',name);
end

function values=expandJoint(values,nj)
values=values(:);
if isscalar(values), values=repmat(values,nj,1); end
assert(numel(values)==nj,'Expected scalar or one value per joint.');
end

function ratio=relative(numerator,denominator)
ratio=NaN;
if isfinite(denominator)&&abs(denominator)>1e-12, ratio=numerator/denominator; end
end

function out=appendTable(out,rows)
if isempty(out), out=rows; else, out=[out;rows]; end
end

function makeOverview(series,summary,annual,files)
colors=[.15 .45 .75;.90 .50 .10;.30 .65 .40;.60 .30 .70];
labels={'K_{sat}=10^{-7}, fixed WT','K_{sat}=10^{-9}, fixed WT', ...
 'K_{sat}=10^{-10}, fixed WT','K_{sat}=10^{-10}, seasonal WT'};
fig=figure('Visible','off','Color','w','Position',[80 80 1400 900]);
cleanFig=onCleanup(@() close(fig)); %#ok<NASGU>
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); hold(ax,'on');
for k=1:4
 d=series{k}.daily;
 plot(ax,d.day_of_year_binEnd-.5,d.actual_L_per_day_per_m, ...
  'Color',colors(k,:),'LineWidth',1.7,'DisplayName',labels{k});
end
title(ax,'(a) Year-20 daily actual joint leakage');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'Leakage (L day^{-1} m^{-1})');
xlim(ax,[0 365]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);

ax=nexttile(t);
matrix=[summary.Q20_actual_L_per_m(1:4),summary.Q20_sameT_static_L_per_m(1:4), ...
 summary.Q20_sameT_seasonal_L_per_m(1:4),summary.Q20_independentT_seasonal_L_per_m(1:4)];
bars=bar(ax,1:4,matrix,'grouped');
barColors=[.1 .35 .65;.75 .75 .75;.90 .55 .20;.25 .60 .35];
for k=1:4, bars(k).FaceColor=barColors(k,:); end
xticks(ax,1:4); xticklabels(ax,{'Base','Transition','Supply','Supply + WT'});
ylabel(ax,'Year-20 annual leakage (L m^{-1})');
title(ax,'(b) Coupled leakage and head approximations');
legend(ax,bars,{'Actual TH','Same T, fixed head','Same T, actual WT', ...
 'Independent T, actual WT'},'Location','best','FontSize',9); formatAxes(ax);

ax=nexttile(t); hold(ax,'on');
semilogx(ax,summary.Ksat_m_per_s(1:3),summary.actual_over_sameT_seasonal(1:3), ...
 '-o','Color',[.15 .45 .75],'LineWidth',1.6,'MarkerSize',7);
yline(ax,1,'k:','No supply loss','LabelHorizontalAlignment','left');
xlabel(ax,'Soil K_{sat} referenced at 20 ^{\circ}C (m s^{-1})');
ylabel(ax,'Q_{actual}/Q_{sameT,referenceHead}');
title(ax,'(c) Computed fixed-water-table permeability scan');
ylim(ax,[min(0,min(summary.actual_over_sameT_seasonal(1:3))-.05),1.05]);
formatAxes(ax); set(ax,'XScale','log');

ax=nexttile(t); hold(ax,'on');
for k=1:4
 d=annual(annual.caseName==string(series{k}.name),:);
 plot(ax,d.year,d.annualActual_L_per_m,'-o','Color',colors(k,:), ...
  'LineWidth',1.4,'MarkerSize',3,'DisplayName',labels{k});
end
title(ax,'(d) Complete transient annual leakage history');
xlabel(ax,'Simulation year'); ylabel(ax,'Annual leakage (L m^{-1})');
xlim(ax,[1 20]); xticks(ax,[1 5 10 15 20]);
legend(ax,'Location','best','FontSize',9); formatAxes(ax);
title(t,'Illustrative freeze-thaw joint leakage: soil supply versus joint resistance', ...
 'FontName','Arial','FontSize',15,'FontWeight','normal');
exportgraphics(fig,files.overviewPNG,'Resolution',300);
exportgraphics(fig,files.overviewPDF,'ContentType','vector');
end

function makeSeasonal(actual,background,files)
a=actual.daily; b=background.daily;
fig=figure('Visible','off','Color','w','Position',[80 80 1400 900]);
cleanFig=onCleanup(@() close(fig)); %#ok<NASGU>
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); hold(ax,'on');
plot(ax,a.day_of_year_binEnd-.5,a.actual_L_per_day_per_m, ...
 'Color',[.10 .35 .65],'LineWidth',1.8,'DisplayName','Actual coupled TH');
plot(ax,a.day_of_year_binEnd-.5,a.sameT_seasonal_L_per_day_per_m, ...
 '--','Color',[.90 .50 .10],'LineWidth',1.5,'DisplayName','Same T, reference head');
plot(ax,b.day_of_year_binEnd-.5,b.sameT_seasonal_L_per_day_per_m, ...
 ':','Color',[.25 .60 .35],'LineWidth',1.8,'DisplayName','Independent T, reference head');
title(ax,'(a) Matched seasonal-water-table leakage comparison');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'Leakage (L day^{-1} m^{-1})');
xlim(ax,[0 365]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);

ax=nexttile(t); hold(ax,'on');
plot(ax,a.day_of_year_binEnd-.5,a.seasonalReferenceHead_mean_m,'--', ...
 'Color',[.90 .50 .10],'LineWidth',1.7,'DisplayName','Reference head, actual WT');
plot(ax,a.day_of_year_binEnd-.5,a.actualDrivingHead_mean_m, ...
 'Color',[.10 .35 .65],'LineWidth',1.7,'DisplayName','Actual joint driving head');
plot(ax,a.day_of_year_binEnd-.5,a.staticReferenceHead_mean_m,':', ...
 'Color',[.5 .5 .5],'LineWidth',1.3,'DisplayName','Fixed-WT reference head');
title(ax,'(b) Supply-induced head loss at the active joint');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'Driving head (m; before entry threshold)');
xlim(ax,[0 365]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);

ax=nexttile(t); hold(ax,'on');
plot(ax,a.day_of_year_binEnd-.5,a.factor_mean,'Color',[.10 .35 .65], ...
 'LineWidth',1.8,'DisplayName','Actual coupled TH');
plot(ax,b.day_of_year_binEnd-.5,b.factor_mean,':','Color',[.25 .60 .35], ...
 'LineWidth',1.8,'DisplayName','Independent no-extraction TH');
title(ax,'(c) Thermal inputs compared independently');
xlabel(ax,'Day of repeated 365-day year'); ylabel(ax,'C_{eff}/C_{ref} (time mean)');
xlim(ax,[0 365]); legend(ax,'Location','best','FontSize',9); formatAxes(ax);

ax=nexttile(t);
v=[actual.monthly.actual_L_per_m,actual.monthly.sameT_seasonal_L_per_m, ...
 background.monthly.sameT_seasonal_L_per_m];
bars=bar(ax,1:12,v,'grouped');
barColors=[.10 .35 .65;.90 .50 .10;.25 .60 .35];
for k=1:3, bars(k).FaceColor=barColors(k,:); end
xticks(ax,1:12); xticklabels(ax,{'Jan','Feb','Mar','Apr','May','Jun', ...
 'Jul','Aug','Sep','Oct','Nov','Dec'});
title(ax,'(d) Conservative calendar-month integrals');
xlabel(ax,'Calendar month'); ylabel(ax,'Monthly leakage (L m^{-1})');
legend(ax,bars,{'Actual TH','Same T, reference head','Independent T, reference head'}, ...
 'Location','best','FontSize',9); formatAxes(ax);
title(t,'Low-permeability illustrative case: seasonal head, icing and soil supply', ...
 'FontName','Arial','FontSize',15,'FontWeight','normal');
exportgraphics(fig,files.seasonalPNG,'Resolution',300);
exportgraphics(fig,files.seasonalPDF,'ContentType','vector');
end

function formatAxes(ax)
set(ax,'FontName','Arial','FontSize',11,'LineWidth',.8,'Box','on','TickDir','out');
grid(ax,'on'); ax.GridAlpha=.15;
end
