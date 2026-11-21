function fig=plot_joint_only_year20(result,outputFile)
%PLOT_JOINT_ONLY_YEAR20 Daily total and monthly rates by longitudinal joint.
% Volumes are L per metre of tunnel axis, not whole-tunnel volumes.
if nargin<2, outputFile=''; end
if isfield(result,'lastYear'), out=result.lastYear; else, out=result; end
assert(out.parameters.liningK==0 && ...
    out.leakage.bodyVolume_m3_per_m==0 && ...
    abs(out.timeDay/365-20)<1e-8, ...
    'Supply a completed age-20 joint-only result.');
L=out.leakage;
fig=figure('Name','Year-20 longitudinal joint leakage', ...
    'Color','w','Visible','off','Position',[100 100 1100 680]);
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(1:365,1000*L.dailyTotalVolume_m3_per_m,'k-','LineWidth',1.3);
xlim([1 365]); grid on;
xlabel('Day of year 20'); ylabel('Leakage (L/day/m)');
title(sprintf('Year 20: %.2f L/m through longitudinal joints', ...
    1000*L.totalVolume_m3_per_m));
nexttile;
active=find(any(L.monthlyVolumeByJoint_m3_per_m>0,1));
if isempty(active)
    bar(1:12,zeros(12,1));
else
    bar(1:12,1000*L.monthlyVolumeByJoint_m3_per_m(:,active),'stacked');
end
xlim([.5 12.5]); grid on;
xticks(1:12); xlabel('Month'); ylabel('Leakage (L/month/m)');
if ~isempty(active)
    legend(compose('Joint %d',active),'Location','eastoutside');
end
if ~isempty(outputFile)
    exportgraphics(fig,outputFile,'Resolution',200);
end
end
