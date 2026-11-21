function [Sfine,report]=refine_active_joint_mesh(S,angles,widths,level)
%REFINE_ACTIVE_JOINT_MESH Resolve selected longitudinal joints on both sides.
% [SFINE,REPORT]=REFINE_ACTIVE_JOINT_MESH(S,ANGLES,WIDTHS,LEVEL) inserts
% circular lining nodes and conforming soil interface nodes at the centre
% and ends of each drainage band. The band arc spacing is WIDTH/2^LEVEL.
% All radial layers of the circular Q4 shell use the same angular cuts.
% A graded soil patch is retriangulated with constrained Delaunay, retaining
% original soil points except points swallowed by the refined hole or too
% close to a hole edge to give a nonnegative boundary cotangent weight.
%
% This routine returns structs only and never writes a mesh or a field.
% New circular vertices replace coarse chords, so the local polygonal
% domain changes slightly. Start a fresh solution on SFINE; an old state
% cannot simply be copied to it. Angles are ATAN2 about [S.xc,S.yc].
% Refinement of shell sectors attached to invert elements is deliberately
% rejected, since the attached elements would otherwise acquire hanging
% nodes. The default leaking crown joint is supported by the supplied mesh.
required={'pAll','tSoil','qTunnelGlobal','soilTunnelCoupleTable', ...
    'tunnelInnerRingGlobal','xc','yc'};
assert(isstruct(S)&&isscalar(S)&&all(isfield(S,required)), ...
    'A circular shell mesh with a soil-lining coupling table is required.');
assert(isnumeric(angles)&&isreal(angles)&&isvector(angles)&& ...
    ~isempty(angles)&&all(isfinite(angles(:))),'Supply active joint angles.');
assert(isnumeric(widths)&&isreal(widths)&& ...
    all(isfinite(widths(:)))&&all(widths(:)>0),'Widths must be positive.');
assert(isscalar(level)&&isfinite(level)&&level==fix(level)&& ...
    level>=1&&level<=5,'Level must be an integer from 1 to 5.');
angles=mod(double(angles(:)),2*pi); widths=double(widths(:));
if isscalar(widths), widths=repmat(widths,numel(angles),1); end
assert(numel(widths)==numel(angles),'Supply one width per active joint.');
xy=double(S.pAll); xyOriginal=xy; Nold=size(xy,1);
assert(size(xy,2)==2&&all(isfinite(xy(:))),'Invalid mesh coordinates.');
triOld=double(S.tSoil); qOld=double(S.qTunnelGlobal);
validateConnectivity(triOld,3,Nold,'tSoil');
validateConnectivity(qOld,4,Nold,'qTunnelGlobal');
center=double([S.xc S.yc]); tol=1e-10; angleTol=1e-10;
shellIDs=unique(qOld(:)); delta=xy(shellIDs,:)-center;
shellR=hypot(delta(:,1),delta(:,2)); shellA=mod(atan2(delta(:,2),delta(:,1)),2*pi);
shellA(abs(shellA-2*pi)<angleTol)=0;
radii=sort(unique(round(shellR*1e8)/1e8),'descend');
assert(numel(radii)>=2&&all(radii>0),'Shell must have circular radial layers.');
outerRadius=radii(1);
assert(all(widths<pi*outerRadius/4),'Joint refinement bands are too wide.');
for i=1:numel(angles)
    for j=i+1:numel(angles)
        distance=abs(wrapAngle(angles(i)-angles(j)))*outerRadius;
        assert(distance>(widths(i)+widths(j))/2+tol, ...
            'Active joint bands overlap.');
    end
end
originalAngles=shellA(abs(shellR-outerRadius)<1e-7);
originalAngles=sort(originalAngles);
nold=numel(originalAngles); nr=numel(radii);
assert(nold>=8&&numel(shellIDs)==nr*nold&& ...
    size(qOld,1)==(nr-1)*nold, ...
    'Only a complete conforming circular annular Q4 shell is supported.');
ringOriginal=zeros(nr,nold);
for r=1:nr
    ids=shellIDs(abs(shellR-radii(r))<1e-7);
    a=shellA(abs(shellR-radii(r))<1e-7);
    assert(numel(ids)==nold,'Inconsistent shell ring node count.');
    for k=1:nold
        match=find(abs(wrapAngle(a-originalAngles(k)))<angleTol);
        assert(isscalar(match),'Shell rings do not share angular cuts.');
        ringOriginal(r,k)=ids(match);
    end
end
assert(isequal(sort(unique(S.tunnelInnerRingGlobal(:))), ...
    sort(ringOriginal(end,:)')),'The inner ring does not match the shell.');

% Each selected band has exact end cuts; outside it, intervals grow by two
% until the original angular vertices enclosing that band are reached.
newAngles=originalAngles;
for j=1:numel(angles)
    relative=sort(wrapAngle(originalAngles-angles(j))*outerRadius);
    half=widths(j)/2; dx=widths(j)/2^level;
    left=max(relative(relative<=-half+tol));
    right=min(relative(relative>=half-tol));
    assert(~isempty(left)&&~isempty(right),'Cannot bracket a selected band.');
    cuts=linspace(-half,half,2^level+1)';
    cuts=[cuts;gradedCuts(-half,left,dx,-1);gradedCuts(half,right,dx,1)]; %#ok<AGROW>
    newAngles=[newAngles;mod(angles(j)+cuts/outerRadius,2*pi)]; %#ok<AGROW>
end
newAngles=angularUnique(newAngles,angleTol); nnew=numel(newAngles);
isOld=false(nnew,1); oldAngleIndex=zeros(nnew,1);
for k=1:nnew
    match=find(abs(wrapAngle(originalAngles-newAngles(k)))<angleTol);
    if ~isempty(match)
        assert(isscalar(match)); isOld(k)=true; oldAngleIndex(k)=match;
    end
end

% Refining an angular sector must not split a shared shell/invert edge.
attached=zeros(0,1);
for name={'qInvertGlobal','tInvertTriGlobal'}
    if isfield(S,name{1})&&~isempty(S.(name{1}))
        attached=[attached;S.(name{1})(:)]; %#ok<AGROW>
    end
end
attached=unique(attached);
for k=find(~isOld)'
    sector=find(originalAngles<newAngles(k)-angleTol,1,'last');
    if isempty(sector), sector=nold; end
    next=mod(sector,nold)+1;
    endpoints=ringOriginal(:,[sector next]);
    assert(isempty(intersect(endpoints(:),attached)), ...
        ['Selected sector shares nodes with invert concrete; conforming ' ...
         'invert refinement is not implemented. Select a crown joint.']);
end

ring=zeros(nr,nnew);
for r=1:nr
    for k=1:nnew
        if isOld(k)
            ring(r,k)=ringOriginal(r,oldAngleIndex(k));
        else
            xy(end+1,:)=center+radii(r)*[cos(newAngles(k)),sin(newAngles(k))]; %#ok<AGROW>
            ring(r,k)=size(xy,1);
        end
    end
end
qNew=zeros((nr-1)*nnew,4);
for r=1:nr-1
    next=[2:nnew 1]; rows=(r-1)*nnew+(1:nnew);
    qNew(rows,:)=[ring(r,:)',ring(r,next)',ring(r+1,next)',ring(r+1,:)'];
end

% Soil and concrete interface coordinates intentionally have distinct IDs.
soilIDs=unique(triOld(:)); soilBoundary=boundaryEdges(triOld);
interfaceOld=sort(double(S.soilTunnelCoupleTable(:,2:3)),2);
assert(size(unique(interfaceOld,'rows'),1)==nold, ...
    'The soil-lining interface does not match the circular shell.');
assert(all(ismember(interfaceOld,soilBoundary,'rows')), ...
    'Original interface edges are not soil boundary edges.');
outerBoundary=soilBoundary(~ismember(soilBoundary,interfaceOld,'rows'),:);
soilOuter=zeros(nnew,1);
for k=1:nnew
    point=xy(ring(1,k),:);
    [xy,soilIDs,id]=addSoilPoint(xy,soilIDs,point,tol);
    soilOuter(k)=id;
end
holeEdges=[soilOuter,soilOuter([2:nnew 1])];
holePolygon=xy(soilOuter,:);
oldLocal=zeros(Nold,1); oldLocal(unique(triOld(:)))=1:numel(unique(triOld(:)));
originalTR=triangulation(oldLocal(triOld),xyOriginal(unique(triOld(:)),:));

% Near each active band create radial layers with comparable radial/arc
% spacing. Coarsen angular spacing away from the interface; keep the patch
% local to a 1 m radial buffer and its original enclosing angular vertices.
for j=1:numel(angles)
    dx=widths(j)/2^level; distance=dx; lastDistance=0; nlayer=0;
    while distance<=1+tol
        nlayer=nlayer+1; radialStep=distance-lastDistance;
        spacing=max(dx,radialStep); halfPatch=widths(j)/2+max(.1,2*distance);
        offset=(-ceil(halfPatch/spacing):ceil(halfPatch/spacing))'*spacing;
        offset=offset(abs(offset)<=halfPatch+tol);
        localA=angles(j)+offset/(outerRadius+distance);
        points=center+(outerRadius+distance)*[cos(localA),sin(localA)];
        inOriginal=~isnan(pointLocation(originalTR,points));
        for k=find(inOriginal)'
            [xy,soilIDs]=addSoilPoint(xy,soilIDs,points(k,:),tol);
        end
        lastDistance=distance;
        distance=distance+dx*2^nlayer;
    end
end
% Midpoint guards keep long unchanged hole edges from acquiring obtuse
% boundary triangles during global retriangulation. They are only added
% where the original soil domain contains them.
a=xy(holeEdges(:,1),:); b=xy(holeEdges(:,2),:);
edgeLength=hypot(b(:,1)-a(:,1),b(:,2)-a(:,2)); midpoint=(a+b)/2;
normal=midpoint-center; normal=normal./hypot(normal(:,1),normal(:,2));
guards=midpoint+.65*edgeLength.*normal;
inOriginal=~isnan(pointLocation(originalTR,guards));
for k=find(inOriginal)'
    [xy,soilIDs]=addSoilPoint(xy,soilIDs,guards(k,:),tol);
end

% Remove soil points inside the expanded polygon and in the exterior half
% of a hole edge's diameter disk. No interface or outer boundary is removed.
protected=unique([soilOuter;outerBoundary(:)]);
[inside,on]=inpolygon(xy(soilIDs,1),xy(soilIDs,2),holePolygon(:,1),holePolygon(:,2));
remove=inside&~on&~ismember(soilIDs,protected);
for k=1:size(holeEdges,1)
    d2=sum((xy(soilIDs,:)-midpoint(k,:)).^2,2);
    remove=remove|((d2<(edgeLength(k)/2)^2-1e-12) ...
        &~ismember(soilIDs,protected));
end
removedIDs=soilIDs(remove); soilIDs=soilIDs(~remove);
local=zeros(size(xy,1),1); local(soilIDs)=1:numel(soilIDs);
constraints=unique(sort([outerBoundary;holeEdges],2),'rows');
assert(all(local(constraints(:))>0),'A constrained boundary point was removed.');
DT=delaunayTriangulation(xy(soilIDs,:),local(constraints));
assert(size(DT.Points,1)==numel(soilIDs), ...
    'Delaunay merged duplicate soil points unexpectedly.');
tri=soilIDs(DT.ConnectivityList);
centroid=(xy(tri(:,1),:)+xy(tri(:,2),:)+xy(tri(:,3),:))/3;
[inHole,onHole]=inpolygon(centroid(:,1),centroid(:,2),holePolygon(:,1),holePolygon(:,2));
inOriginal=~isnan(pointLocation(originalTR,centroid));
tri=tri(inOriginal&~inHole&~onHole,:);
signed=signedArea2(xy,tri);
assert(all(abs(signed)>1e-12),'A degenerate refined soil triangle was created.');
tri(signed<0,[2 3])=tri(signed<0,[3 2]);
newBoundary=boundaryEdges(tri);
assert(isequal(newBoundary,constraints), ...
    'Refined soil has a missing constraint, gap, or extra boundary.');
assert(size(unique(sort(tri,2),'rows'),1)==size(tri,1), ...
    'Refined soil contains duplicate triangles.');

% Rebuild the coupling table, including containing triangle and local edge.
couple=zeros(nnew,9);
for k=1:nnew
    ends=holeEdges(k,:);
    containing=find(sum(ismember(tri,ends),2)==2);
    assert(isscalar(containing),'Interface edge must have one soil triangle.');
    pairs=[tri(containing,[1 2]);tri(containing,[2 3]);tri(containing,[3 1])];
    edge=find(all(sort(pairs,2)==sort(ends),2));
    assert(isscalar(edge),'Cannot identify the soil triangle local edge.');
    couple(k,:)=[k,ends,containing,edge,ring(1,k),ring(1,mod(k,nnew)+1),k,1];
end
Sfine=S; Sfine.pAll=xy; Sfine.tSoil=tri;
Sfine.qTunnelGlobal=qNew; Sfine.tunnelInnerRingGlobal=ring(end,:)';
Sfine.soilTunnelCoupleTable=couple;

% Remove all unused points, retaining the distinction between soil and
% coincident lining nodes. Remap every structural/boundary connectivity.
fields={'qTunnelGlobal','qInvertGlobal','tInvertTriGlobal', ...
    'tunnelInnerRingGlobal','invertBottomRingGlobal', ...
    'invertTopRingGlobal','bInvertTopGlobal'};
used=unique(tri(:));
for name=fields
    if isfield(Sfine,name{1}), used=union(used,Sfine.(name{1})(:)); end
end
map=zeros(size(xy,1),1); map(used)=1:numel(used);
Sfine.pAll=xy(used,:); Sfine.tSoil=map(tri);
for name=fields
    if isfield(Sfine,name{1})&&~isempty(Sfine.(name{1}))
        value=Sfine.(name{1}); assert(all(map(value(:))>0));
        Sfine.(name{1})=reshape(map(value),size(value));
    end
end
for column=[2 3 6 7], couple(:,column)=map(couple(:,column)); end
Sfine.soilTunnelCoupleTable=couple;
originalNode=zeros(numel(used),1); originalNode(used<=Nold)=used(used<=Nold);
Sfine.newNodeToOriginalNode=originalNode;

actualLength=zeros(numel(angles),2); count=zeros(numel(angles),1);
for j=1:numel(angles)
    relative=wrapAngle(newAngles-angles(j))*outerRadius;
    nextRelative=relative+wrapAngle(newAngles([2:nnew 1])-newAngles)*outerRadius;
    low=min(relative,nextRelative); high=max(relative,nextRelative);
    overlap=high>-widths(j)/2+tol&low<widths(j)/2-tol;
    lengths=edgeLength(overlap); assert(~isempty(lengths));
    actualLength(j,:)=[min(lengths),max(lengths)]; count(j)=numel(lengths);
    assert(max(lengths)<=widths(j)/2^level*(1+1e-6), ...
        'Selected drainage band is not resolved at the requested spacing.');
end
oldArea=sum(abs(signedArea2(xyOriginal,triOld)))/2;
newArea=sum(abs(signedArea2(Sfine.pAll,Sfine.tSoil)))/2;
report=struct('level',level,'activeAnglesRad',angles,'jointWidths_m',widths, ...
    'requestedBandSpacing_m',widths/2^level, ...
    'actualBandEdgeLengthRange_m',actualLength,'bandInterfaceEdges',count, ...
    'originalNodes',Nold,'refinedNodes',size(Sfine.pAll,1), ...
    'originalSoilTriangles',size(triOld,1),'refinedSoilTriangles',size(tri,1), ...
    'originalShellQuads',size(qOld,1),'refinedShellQuads',size(qNew,1), ...
    'originalInterfaceEdges',nold,'refinedInterfaceEdges',nnew, ...
    'originalSoilArea_m2',oldArea,'refinedSoilArea_m2',newArea, ...
    'soilAreaChange_m2',newArea-oldArea,'removedSoilPoints',numel(removedIDs), ...
    'originalMinimumSoilAngle_deg',min(triangleMinAngle(xyOriginal,triOld)), ...
    'refinedMinimumSoilAngle_deg',min(triangleMinAngle(Sfine.pAll,Sfine.tSoil)), ...
    'unusedNodes',size(xy,1)-numel(used),'retainedUnusedNodes',0, ...
    'freshInitialStateRequired',true,'method','graded circular Q4 + constrained Delaunay');
[negativeAll,negativeInterface]=negativeHydraulicEdges(Sfine.pAll, ...
    Sfine.tSoil,couple(:,2:3));
report.negativeHydraulicEdges=negativeAll;
report.negativeInterfaceHydraulicEdges=negativeInterface;
assert(negativeInterface==0,'Refined interface has negative hydraulic edge weights.');
Sfine.jointRefinementReport=report;
end

function validateConnectivity(t,width,N,name)
assert(isnumeric(t)&&isreal(t)&&size(t,2)==width&&all(isfinite(t(:))) ...
    &&all(t(:)==fix(t(:)))&&all(t(:)>=1&t(:)<=N), ...
    'Invalid %s connectivity.',name);
assert(all(sum(diff(sort(t,2),1,2)==0,2)==0),'Repeated element nodes in %s.',name);
end

function cuts=gradedCuts(start,limit,spacing,direction)
cuts=zeros(0,1); position=start;
while direction*(limit-position)>spacing+1e-10
    position=position+direction*spacing; cuts(end+1,1)=position; %#ok<AGROW>
    spacing=2*spacing;
end
end

function a=angularUnique(a,tol)
a=sort(mod(a(:),2*pi)); a(abs(a-2*pi)<tol)=0; a=sort(a);
a=a([true;diff(a)>tol]);
if numel(a)>1&&2*pi-a(end)+a(1)<tol, a(end)=[]; end
end

function a=wrapAngle(a)
a=atan2(sin(a),cos(a));
end

function [xy,soilIDs,id]=addSoilPoint(xy,soilIDs,point,tol)
distance2=sum((xy(soilIDs,:)-point).^2,2);
[d,index]=min(distance2);
if d<=tol^2
    id=soilIDs(index);
else
    xy(end+1,:)=point; id=size(xy,1); soilIDs(end+1,1)=id;
end
end

function b=boundaryEdges(t)
e=sort([t(:,[1 2]);t(:,[2 3]);t(:,[3 1])],2);
[u,~,g]=unique(e,'rows'); count=accumarray(g,1);
assert(all(count<=2),'Nonmanifold soil edge.'); b=u(count==1,:);
end

function a=signedArea2(x,t)
p=x(t(:,1),:); q=x(t(:,2),:); r=x(t(:,3),:);
a=(q(:,1)-p(:,1)).*(r(:,2)-p(:,2))-(q(:,2)-p(:,2)).*(r(:,1)-p(:,1));
end

function angle=triangleMinAngle(x,t)
p=x(t(:,1),:); q=x(t(:,2),:); r=x(t(:,3),:);
a=abs(signedArea2(x,t));
l12=hypot(p(:,1)-q(:,1),p(:,2)-q(:,2));
l23=hypot(q(:,1)-r(:,1),q(:,2)-r(:,2));
l31=hypot(r(:,1)-p(:,1),r(:,2)-p(:,2));
cosines=[(l12.^2+l31.^2-l23.^2)./(2*l12.*l31), ...
    (l12.^2+l23.^2-l31.^2)./(2*l12.*l23), ...
    (l23.^2+l31.^2-l12.^2)./(2*l23.*l31)];
assert(all(a>1e-12)); angle=min(acosd(max(-1,min(1,cosines))),[],2);
end

function [negativeAll,negativeInterface]=negativeHydraulicEdges(x,t,interface)
p=x(t(:,1),:); q=x(t(:,2),:); r=x(t(:,3),:); a2=signedArea2(x,t);
g1=[q(:,2)-r(:,2),r(:,1)-q(:,1)]./a2;
g2=[r(:,2)-p(:,2),p(:,1)-r(:,1)]./a2;
g3=[p(:,2)-q(:,2),q(:,1)-p(:,1)]./a2;
area=abs(a2)/2;
w=[-area.*sum(g1.*g2,2);-area.*sum(g2.*g3,2);-area.*sum(g3.*g1,2)];
e=sort([t(:,[1 2]);t(:,[2 3]);t(:,[3 1])],2);
[u,~,g]=unique(e,'rows'); w=accumarray(g,w);
[found,index]=ismember(sort(interface,2),u,'rows'); assert(all(found));
negativeAll=nnz(w< -1e-12); negativeInterface=nnz(w(index)< -1e-12);
end
