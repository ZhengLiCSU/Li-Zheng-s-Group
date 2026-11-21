function artifacts=plot_thm20_time_step_verification(caseFolder)
%PLOT_THM20_TIME_STEP_VERIFICATION Plot saved target-year checks, no solving.
% The comparison is conditional on the same year-19 state. Monthly errors
% use the fine ANNUAL volume as denominator, as in the pre-set criterion.
if nargin<1||isempty(caseFolder)
    caseFolder=fullfile(fileparts(mfilename('fullpath')),'thm20_20260930_analytic_p1');
end
outputFolder=fullfile(caseFolder,'year20_time_step_verification');
saved=load(fullfile(outputFolder,'year20_time_step_report.mat'),'report');
r=saved.report; m=r.monthly;
assert(height(m)==12&&isequal(m.month,(1:12)'), ...
    'THM20Plot:Months','The saved comparison must contain all twelve months.');
assert(all(isfinite(m.reference_L_per_m))&&all(isfinite(m.fine_L_per_m)) ...
    &&all(m.reference_L_per_m>=0)&&all(m.fine_L_per_m>=0), ...
    'THM20Plot:Volumes','Saved monthly volumes must be finite and nonnegative.');
f=figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[2 2 18 17]);
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
layout=tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');
monthNames={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
blue=[0.13 0.40 0.68]; orange=[0.87 0.47 0.15];

ax=nexttile(layout);
b=bar(ax,m.month,[m.reference_L_per_m m.fine_L_per_m],0.82);
b(1).FaceColor=blue; b(2).FaceColor=orange;
b(1).EdgeColor='none'; b(2).EdgeColor='none';
ylabel(ax,'Monthly leakage (L/m)');
title(ax,{'(a) Computed monthly leakage',sprintf( ...
    'Annual: %.5f vs %.5f L/m; difference / fine annual = %.4f%%', ...
    r.referenceAnnualLeakage_L_per_m,r.fineAnnualLeakage_L_per_m, ...
    r.annualRelativeDifference_percent)},'FontWeight','normal');
legend(ax,{sprintf('Maximum step %.1f day',r.referenceDtMax_day), ...
    sprintf('Maximum step %.1f day',r.fineDtMax_day)}, ...
    'Location','northwest','Box','off');
formatAxis(ax,monthNames);

ax=nexttile(layout);
b=bar(ax,m.month,m.differenceRelativeToFineAnnual_percent,0.6, ...
    'FaceColor',blue,'EdgeColor','none');
threshold=yline(ax,r.leakageAcceptanceLimit_percent,'--', ...
    'Color',orange,'LineWidth',1.5);
ylabel(ax,{'Absolute monthly difference','/ fine annual leakage (%)'});
xlabel(ax,'Calendar month');
title(ax,'(b) Monthly error normalized by the fine annual volume', ...
    'FontWeight','normal');
legend(ax,[b threshold],{'Monthly difference',sprintf( ...
    'Pre-set %.0f%% criterion',r.leakageAcceptanceLimit_percent)}, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
ylim(ax,[0 1.18*max([r.leakageAcceptanceLimit_percent; ...
    m.differenceRelativeToFineAnnual_percent])]);
formatAxis(ax,monthNames);
title(layout,'Target-year check from the same year-19 THM state', ...
    'FontName','Arial','FontSize',12,'FontWeight','normal');

artifacts=struct('png',fullfile(outputFolder,'year20_time_step_check.png'), ...
    'pdf',fullfile(outputFolder,'year20_time_step_check.pdf'));
exportgraphics(f,artifacts.png,'Resolution',300);
exportgraphics(f,artifacts.pdf,'ContentType','vector');
end

function formatAxis(ax,monthNames)
set(ax,'FontName','Arial','FontSize',10,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XTick',1:12,'XTickLabel',monthNames, ...
    'XLim',[0.4 12.6]);
grid(ax,'on'); ax.GridAlpha=0.13;
end
