function diagnostic=diagnose_leakage_timing(result,outputFile)
%DIAGNOSE_LEAKAGE_TIMING Compare imposed air, lining state, and daily leakage.
% Accepts a complete 365-day v8/v9 solver output or annual-controller result.
if nargin<2, outputFile=''; end
if isfield(result,'lastYear'), out=result.lastYear; else, out=result; end
assert(numel(out.leakage.dailyTotalVolume_m3_per_m)==365 && ...
    numel(out.snapshots)==12,'Supply a complete annual result.');
assert(isfield(out.mesh,'soilTunnelCoupleTable') && ...
    ~isempty(out.mesh.soilTunnelCoupleTable),'The mesh has no lining interface.');
interfaceNodes=unique(out.mesh.soilTunnelCoupleTable(:,[2 3]));
[found,soilIndex]=ismember(interfaceNodes,out.soilNodes);
assert(all(found),'Lining interface nodes are missing from the soil domain.');
p=out.parameters;
if isfield(out.mesh,'xc') && isfield(out.mesh,'yc')
    center=[out.mesh.xc,out.mesh.yc];
else
    center=(min(out.mesh.pAll,[],1)+max(out.mesh.pAll,[],1))/2;
end
outerRadius=hypot(out.mesh.pAll(interfaceNodes,1)-center(1), ...
    out.mesh.pAll(interfaceNodes,2)-center(2));
liningZDrop=zeros(size(outerRadius));
if isfield(out,'version') && startsWith(out.version,'TH_v9')
    liningZDrop=-p.liningThickness*(out.mesh.pAll(interfaceNodes,2)-center(2))./outerRadius;
end
dayMid=(.5:1:364.5)';
phase=2*pi*dayMid/p.yearDays;
topAir=(p.topMax+p.topMin)/2-(p.topMax-p.topMin)/2*cos(phase);
tunnelAir=(p.tunnelMax+p.tunnelMin)/2-(p.tunnelMax-p.tunnelMin)/2*cos(phase);
if ~isempty(p.forcing.topAirC)
    topAir=interp1(p.forcing.day,p.forcing.topAirC,dayMid,'linear');
end
if ~isempty(p.forcing.tunnelAirC)
    tunnelAir=interp1(p.forcing.day,p.forcing.tunnelAirC,dayMid,'linear');
end
monthDays=[31 28 31 30 31 30 31 31 30 31 30 31];
monthEnds=cumsum(monthDays);
monthStarts=[0 monthEnds(1:end-1)];
topMonthly=zeros(12,1); tunnelMonthly=zeros(12,1);
for month=1:12
    ix=(monthStarts(month)+1):monthEnds(month);
    topMonthly(month)=mean(topAir(ix));
    tunnelMonthly(month)=mean(tunnelAir(ix));
end
diagnostic.monthly=table((1:12)',topMonthly,tunnelMonthly, ...
    1000*out.leakage.monthlyTotalVolume_m3_per_m(:), ...
    'VariableNames',{'Month','TopAirMean_C','TunnelAirMean_C','Leak_L_per_m'});
snap=out.snapshots;
snapDay=zeros(12,1); meanInterfaceT=zeros(12,1);
maxInterfaceHead=zeros(12,1); meanInterfaceIce=zeros(12,1);
maxDrivingHead=zeros(12,1);
meanInterfaceLiquid=zeros(12,1);
for k=1:12
    snapDay(k)=snap(k).timeDay-(out.timeDay-365);
    meanInterfaceT(k)=mean(snap(k).T(interfaceNodes));
    maxInterfaceHead(k)=max(snap(k).h(soilIndex));
    maxDrivingHead(k)=max(snap(k).h(soilIndex)-p.leakage.interiorHead+liningZDrop);
    meanInterfaceIce(k)=mean(snap(k).thetaI(soilIndex));
    meanInterfaceLiquid(k)=mean(snap(k).thetaW(soilIndex));
end
diagnostic.interfaceSnapshots=table(snapDay,meanInterfaceT, ...
    maxInterfaceHead,maxDrivingHead,meanInterfaceIce,meanInterfaceLiquid, ...
    'VariableNames',{'Day','MeanLiningSoilT_C','MaxHead_m', ...
    'MaxLeakageDrivingHead_m','MeanIceVolumeFraction','MeanLiquidVolumeFraction'});
diagnostic.dailyLeak_L_per_day_per_m=1000*out.leakage.dailyTotalVolume_m3_per_m(:);
diagnostic.topAir_C=topAir;
diagnostic.tunnelAir_C=tunnelAir;
fig=figure('Name','Seasonal leakage timing diagnostic','Color','w');
tiledlayout(fig,3,1,'TileSpacing','compact');
nexttile;
plot(dayMid,topAir,'Color',[.35 .55 .75]); hold on;
plot(dayMid,tunnelAir,'Color',[.1 .3 .8]);
plot(snapDay,meanInterfaceT,'o-','Color',[.8 .2 .1],'LineWidth',1.2);
yline(0,':k'); grid on; xlim([0 365]);
ylabel('Temperature (C)');
legend('Ground air','Tunnel air','Mean soil at lining','Location','eastoutside');
nexttile;
plot(snapDay,maxDrivingHead,'o-','Color',[.4 .15 .7],'LineWidth',1.2);
yline(0,':k'); grid on; xlim([0 365]);
ylabel('Max leakage driving head (m)');
nexttile;
semilogy(1:365,max(diagnostic.dailyLeak_L_per_day_per_m,1e-10), ...
    'k-','LineWidth',1.1); grid on; xlim([0 365]);
ylabel('Leak (L/day/m)'); xlabel('Day of 365-day year');
for ax=findall(fig,'Type','axes')'
    xticks(ax,(monthStarts+monthEnds)/2);
    xticklabels(ax,{'Jan','Feb','Mar','Apr','May','Jun', ...
        'Jul','Aug','Sep','Oct','Nov','Dec'});
end
diagnostic.figure=fig;
if ~isempty(outputFile), exportgraphics(fig,outputFile,'Resolution',200); end
end
