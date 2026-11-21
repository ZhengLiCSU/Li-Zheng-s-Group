function fig=plot_segment_leakage(result,outputFile)
%PLOT_SEGMENT_LEAKAGE Daily total and calendar-month segment leakage.
% Values are L per metre of tunnel axis, for the final complete 365-day year.
if nargin<2, outputFile=''; end
if isfield(result,'lastYear'), out=result.lastYear; else, out=result; end
assert(isfield(out,'leakage')&&isfield(out.leakage,'monthlyVolumeBySegment_m3_per_m'), ...
    'Supply a complete 365-day segmented TH result.');
L=out.leakage;
fig=figure('Name','Seasonal tunnel lining leakage','Color','w');
tiledlayout(fig,2,1,'TileSpacing','compact');
nexttile;
plot(1:365,1000*L.dailyTotalVolume_m3_per_m,'k-','LineWidth',1.2);
xlim([1 365]); grid on;
xlabel('Day in final year'); ylabel('Total leakage (L/day/m)');
title(sprintf('Year %d: %.3f L/m total',round(out.timeDay/365), ...
    1000*L.totalVolume_m3_per_m));
nexttile;
bar(1:12,1000*L.monthlyVolumeBySegment_m3_per_m,'stacked');
xlim([.5 12.5]); grid on;
xticks(1:12); xlabel('Calendar month'); ylabel('Segment leakage (L/month/m)');
legend(compose('Segment %d',1:L.segmentCount),'Location','eastoutside');
if ~isempty(outputFile), exportgraphics(fig,outputFile,'Resolution',200); end
end
