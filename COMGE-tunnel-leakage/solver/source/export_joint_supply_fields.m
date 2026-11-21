function files=export_joint_supply_fields(studyFolder,caseName,snapshotBins)
%EXPORT_JOINT_SUPPLY_FIELDS Local FE temperature/ice/head scientific plots.
% Default: base case, reporting bins 3, 4 and 7 (March/April/July bins).
% Use [1 3 4 7] to include the January bin. These model reporting bins are
% equally spaced by 365/12 days, not Gregorian calendar-month endpoints.
% Every panel gives the actual selected simulated phase day.
%
% Coordinates are tangential s and radial-outward y around joint 7. The
% nominal circular interface is y=0 at s=0; actual polygon boundaries are
% drawn from the mesh. Soil T3 and concrete Q4/T3 retain independent nodes.
% Temperature patches use nodal interpolation. Ice color uses thetaI/phi.
% Drawdown contours use source T3 barycentric interpolation, masked by its
% exact triangle inventory; no Delaunay interpolation crosses the lining.
% Dots mark equivalent thermal path samples, not a resolved microchannel.
root=fileparts(mfilename('fullpath')); addpath(root);
if nargin<1||isempty(studyFolder)
    studyFolder=fullfile(root,'joint_supply_study_20260930');
end
if nargin<2||isempty(caseName), caseName='base'; end
if nargin<3||isempty(snapshotBins), snapshotBins=[3 4 7]; end
studyFolder=char(studyFolder); caseName=char(caseName);
snapshotBins=double(snapshotBins(:)');
assert(all(snapshotBins==fix(snapshotBins))&&all(snapshotBins>=1&snapshotBins<=12) ...
    &&numel(unique(snapshotBins))==numel(snapshotBins)&&~isempty(snapshotBins), ...
    'Select distinct reporting-bin indices from one to twelve.');
sourceFile=fullfile(studyFolder,caseName,'TH_segment_spinup_results.mat');
assert(isfile(sourceFile),'Complete source result is missing: %s.',sourceFile);
loaded=load(sourceFile,'result'); result=loaded.result;
assert(result.completedYears==20&&strcmp(result.lastYear.version,'TH_v14_20260930_joint_supply'), ...
    'Local fields require this case''s own complete twenty-year v14 history.');
out=result.lastYear; S=out.mesh; p=out.parameters; snapshots=out.snapshots;
assert(numel(snapshots)==12,'The complete twelve-snapshot twentieth year is required.');
assert(numel(p.leakage.jointAnglesRad)>=7,'Joint 7 is absent.');
N=size(S.pAll,1); soil=out.soilNodes(:);
assert(isequal(soil,unique(S.tSoil(:))),'Source soil-node inventory differs from its mesh.');
theta=p.leakage.jointAnglesRad(7);
radial=[cos(theta);sin(theta)]; tangent=[-sin(theta);cos(theta)];
interface=unique(S.soilTunnelCoupleTable(:,2:3));
R=mean(hypot(S.pAll(interface,1)-S.xc,S.pAll(interface,2)-S.yc));
origin=[S.xc S.yc]+R*radial';
delta=S.pAll-origin; localXY=[delta*tangent,delta*radial];
window=[-.5 .5;-.5 .8];
nearNode=localXY(:,1)>=window(1,1)&localXY(:,1)<=window(1,2) ...
    &localXY(:,2)>=window(2,1)&localXY(:,2)<=window(2,2);
nearSoil=nearNode(soil);
assert(any(nearSoil),'The selected local window contains no soil nodes.');
soilFaces=S.tSoil(faceIntersects(localXY,S.tSoil,window),:);
quads=zeros(0,4); triangles=zeros(0,3);
for name={'qTunnelGlobal','qInvertGlobal'}
    if isfield(S,name{1})&&~isempty(S.(name{1}))
        candidate=S.(name{1});
        quads=[quads;candidate(faceIntersects(localXY,candidate,window),:)]; %#ok<AGROW>
    end
end
if isfield(S,'tInvertTriGlobal')&&~isempty(S.tInvertTriGlobal)
    candidate=S.tInvertTriGlobal;
    triangles=candidate(faceIntersects(localXY,candidate,window),:);
end
outerEdges=S.soilTunnelCoupleTable(:,2:3);
inner=S.tunnelInnerRingGlobal(:); innerEdges=[inner,inner([2:end 1])];
map=zeros(N,1); map(soil)=1:numel(soil);
TR=triangulation(map(S.tSoil),localXY(soil,:));
[gridX,gridY]=meshgrid(linspace(window(1,1),window(1,2),301), ...
    linspace(window(2,1),window(2,2),301));
[containing,barycentric]=pointLocation(TR,[gridX(:),gridY(:)]);
inSoil=~isnan(containing);
sourceTriangle=TR.ConnectivityList(containing(inSoil),:);
allPhase=[snapshots.timeDay]-19*365;
targets=snapshotBins*365/12;
selected=zeros(size(snapshotBins));
temperatureRange=[inf -inf]; drawdownRange=[inf -inf];
for index=1:numel(targets)
    [~,selected(index)]=min(abs(allPhase-targets(index)));
    snap=snapshots(selected(index));
    assert(numel(snap.T)==N&&numel(snap.h)==numel(soil) ...
        &&numel(snap.thetaI)==numel(soil),'Incomplete local-field snapshot.');
    temperatureRange=[min(temperatureRange(1),min(snap.T(nearNode))), ...
        max(temperatureRange(2),max(snap.T(nearNode)))];
    depth=waterDepth(snap.timeDay,p);
    drawdown=S.pAll(soil,2)-min(S.pAll(:,2))-depth-snap.h(:);
    drawdownRange=[min(drawdownRange(1),min(drawdown(nearSoil))), ...
        max(drawdownRange(2),max(drawdown(nearSoil)))];
end
assert(numel(unique(selected))==numel(selected),'Two target bins selected the same snapshot.');
climT=[floor(temperatureRange(1)),ceil(temperatureRange(2))];
if diff(climT)<1, climT=climT+[-.5 .5]; end
contourLevels=[];
if diff(drawdownRange)>1e-8
    contourLevels=linspace(drawdownRange(1),drawdownRange(2),6);
    contourLevels=contourLevels(2:end-1);
end
paths=build_joint_thermal_paths(S,p.leakage.jointAnglesRad(7), ...
    p.leakage.jointFreeze.quadraturePoints);
pathDelta=paths.sampleXY-origin; pathXY=[pathDelta*tangent,pathDelta*radial];
folder=fullfile(studyFolder,'local_field_exports',caseName);
if ~isfolder(folder), mkdir(folder); end
prefix=matlab.lang.makeValidName(caseName);
files=struct('PNG',fullfile(folder,[prefix '_joint7_local_fields.png']), ...
    'PDF',fullfile(folder,[prefix '_joint7_local_fields.pdf']), ...
    'FIG',fullfile(folder,[prefix '_joint7_local_fields.fig']), ...
    'notesTXT',fullfile(folder,[prefix '_field_interpretation.txt']), ...
    'metadataMAT',fullfile(folder,[prefix '_field_metadata.mat']), ...
    'nodeCSV',{cell(numel(selected),1)});
fig=figure('Visible','off','Color','w','Position',[80 80 560*numel(selected) 1000]);
closeFigure=onCleanup(@() close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,2,numel(selected),'TileSpacing','compact','Padding','compact');
title(layout,sprintf('%s: joint 7 local fields in year 20',caseName), ...
    'FontSize',15,'FontWeight','normal','Interpreter','none');
for index=1:numel(selected)
    snap=snapshots(selected(index)); phaseDay=allPhase(selected(index));
    temperature=snap.T(:); iceFraction=snap.thetaI(:)/p.porosity;
    depth=waterDepth(snap.timeDay,p);
    drawdown=S.pAll(soil,2)-min(S.pAll(:,2))-depth-snap.h(:);
    iceGlobal=nan(N,1); iceGlobal(soil)=iceFraction;
    axT=nexttile(layout,index); hold(axT,'on');
    patch(axT,'Faces',soilFaces,'Vertices',localXY,'FaceVertexCData',temperature, ...
        'FaceColor','interp','EdgeColor','none');
    patchConcrete(axT,quads,triangles,localXY,temperature);
    drawBoundary(axT,outerEdges,localXY,window,[.15 .15 .15],.9);
    drawBoundary(axT,innerEdges,localXY,window,[.15 .15 .15],.9);
    plot(axT,pathXY(:,1),pathXY(:,2),'k:','LineWidth',1);
    plot(axT,pathXY(:,1),pathXY(:,2),'ko','MarkerFaceColor','w', ...
        'MarkerSize',3.6,'LineWidth',.8);
    styleAxes(axT,window,index);
    colormap(axT,turbo(256)); clim(axT,climT);
    cb=colorbar(axT); cb.FontSize=11; cb.Label.FontSize=11;
    cb.Label.String='Temperature (°C)';
    title(axT,sprintf('(%c) T: phase day %.2f','a'+index-1,phaseDay), ...
        'FontSize',12,'FontWeight','normal');
    axI=nexttile(layout,numel(selected)+index); hold(axI,'on');
    patch(axI,'Faces',soilFaces,'Vertices',localXY,'FaceVertexCData',iceGlobal, ...
        'FaceColor','interp','EdgeColor','none');
    patchConcrete(axI,quads,triangles,localXY,[]);
    drawBoundary(axI,outerEdges,localXY,window,[.15 .15 .15],.9);
    drawBoundary(axI,innerEdges,localXY,window,[.15 .15 .15],.9);
    gridDrawdown=nan(numel(gridX),1);
    gridDrawdown(inSoil)=sum(barycentric(inSoil,:).*drawdown(sourceTriangle),2);
    if ~isempty(contourLevels)
        [C,handles]=contour(axI,gridX,gridY,reshape(gridDrawdown,size(gridX)), ...
            contourLevels,'LineColor','k','LineWidth',.85);
        clabel(C,handles,'FontSize',10,'Color','k','LabelSpacing',160);
    end
    styleAxes(axI,window,index);
    xlabel(axI,'Tangential distance s (m)','FontSize',11);
    colormap(axI,parula(256)); clim(axI,[0 1]);
    cb=colorbar(axI); cb.FontSize=11; cb.Label.FontSize=11;
    cb.Label.String='Ice fraction \theta_i/\phi'; cb.Label.Interpreter='tex';
    title(axI,sprintf('(%c) Ice; head-drawdown contours (m)', ...
        'a'+numel(selected)+index-1),'FontSize',12,'FontWeight','normal');
    material=repmat("concrete",N,1); material(soil)="soil";
    pressureHead=nan(N,1); pressureHead(soil)=snap.h(:);
    drawdownGlobal=nan(N,1); drawdownGlobal(soil)=drawdown;
    nodeID=(1:N)'; absoluteDay=repmat(snap.timeDay,N,1); phaseDayColumn=repmat(phaseDay,N,1);
    tableData=table(nodeID,material,absoluteDay,phaseDayColumn,S.pAll(:,1),S.pAll(:,2), ...
        localXY(:,1),localXY(:,2),temperature,iceGlobal,pressureHead,drawdownGlobal,nearNode, ...
        'VariableNames',{'globalNodeID','material','absoluteDay','phaseDay', ...
        'globalX_m','globalZ_down_m','tangential_s_m','radialOutward_y_m', ...
        'temperature_C','iceFraction_thetaI_over_porosity','pressureHead_m', ...
        'headDrawdownFromInstantaneousHydrostatic_m','withinFigureWindow'});
    files.nodeCSV{index}=fullfile(folder,sprintf('%s_snapshot%02d_nodes.csv',prefix,selected(index)));
    writetable(tableData,files.nodeCSV{index});
end
exportgraphics(fig,files.PNG,'Resolution',300);
exportgraphics(fig,files.PDF,'ContentType','vector');
savefig(fig,files.FIG);
metadata=struct('sourceFile',sourceFile,'caseName',caseName,'selectedSnapshotIndices',selected, ...
    'actualPhaseDays',allPhase(selected),'requestedReportingBins',snapshotBins, ...
    'originGlobal_m',origin,'radialDirection',radial,'tangentDirection',tangent, ...
    'figureWindow_m',window,'thermalPathSampleGlobalXY_m',paths.sampleXY, ...
    'thermalPathSampleLocalXY_m',pathXY,'temperatureColorLimits_C',climT, ...
    'drawdownContourLevels_m',contourLevels,'waterTablePolicy','actual snapshot forcing', ...
    'figureInterpolation','source T3/Q4 nodal patch interpolation; source T3 barycentric contours');
save(files.metadataMAT,'metadata','-v7');
fid=fopen(files.notesTXT,'w'); assert(fid>=0,'Cannot write field-figure interpretation.');
closeText=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,['Case %s, own complete v14 year-20 fields.\n' ...
    'Actual phase days: %s. Reporting bins are 365/12 days, not calendar month ends.\n' ...
    'Top: independent soil/lining FE nodal temperature interpolation.\n' ...
    'Bottom: stored soil ice content divided by porosity; gray concrete has no soil ice.\n' ...
    'Black lower contours: pressure-head drawdown relative to instantaneous hydrostatic head, in m.\n' ...
    'Contours use barycentric interpolation on existing soil triangles; air/lining gaps are masked.\n' ...
    'Black upper dots: equivalent radial thermal-path sample points, not a resolved water-channel width.\n' ...
    'Local y is outward from the nominal outer circle at joint 7; actual polygon boundaries are drawn.\n' ...
    'CSV preserves full global node IDs, including coincident soil/concrete nodes with separate fields.\n'], ...
    caseName,mat2str(allPhase(selected),7));
fprintf('Local field exports: %s.\n',folder);
end

function keep=faceIntersects(xy,faces,window)
x=xy(:,1); y=xy(:,2);
xf=reshape(x(faces),size(faces)); yf=reshape(y(faces),size(faces));
keep=max(xf,[],2)>=window(1,1)&min(xf,[],2)<=window(1,2) ...
    &max(yf,[],2)>=window(2,1)&min(yf,[],2)<=window(2,2);
end

function patchConcrete(ax,quads,triangles,xy,temperature)
for faces={quads,triangles}
    if isempty(faces{1}), continue; end
    if isempty(temperature)
        patch(ax,'Faces',faces{1},'Vertices',xy,'FaceColor',[.8 .8 .8],'EdgeColor','none');
    else
        patch(ax,'Faces',faces{1},'Vertices',xy,'FaceVertexCData',temperature, ...
            'FaceColor','interp','EdgeColor','none');
    end
end
end

function drawBoundary(ax,edges,xy,window,color,width)
take=faceIntersects(xy,edges,window); edge=edges(take,:);
x=[xy(edge(:,1),1),xy(edge(:,2),1),nan(size(edge,1),1)]';
y=[xy(edge(:,1),2),xy(edge(:,2),2),nan(size(edge,1),1)]';
plot(ax,x(:),y(:),'Color',color,'LineWidth',width);
end

function styleAxes(ax,window,column)
axis(ax,'equal'); xlim(ax,window(1,:)); ylim(ax,window(2,:));
set(ax,'FontSize',11,'LineWidth',.8,'Box','on','YDir','normal','Layer','top','Color','w');
if column==1, ylabel(ax,'Radial outward distance y (m)','FontSize',11); end
end

function depth=waterDepth(day,p)
depth=p.waterTableDepth;
if ~isempty(p.forcing.waterTableDepth_m)
    depth=interp1(p.forcing.day,p.forcing.waterTableDepth_m,mod(day,p.yearDays),'linear');
end
end
