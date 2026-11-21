function report = coarsen_farfield_mesh(sourceFile,outputFile,protectRadius,spacing,minAngleTarget)
%COARSEN_FARFIELD_MESH Trial remeshing outside a protected tunnel zone.
% All original soil triangles touching the protected radial zone or the
% soil-domain boundary, and all concrete/invert elements are retained.
% The remaining soil is retriangulated with constrained Delaunay. This is
% an experimental mesh for verified A/B studies, not a default replacement.
root = fileparts(mfilename('fullpath'));
if nargin<1 || isempty(sourceFile)
    sourceFile=fullfile(root,'example_mesh_repaired_interface.mat');
end
if nargin<2 || isempty(outputFile)
    outputFile=fullfile(root,'farfield_mesh_trial_20260923', ...
        'example_mesh_farfield_coarse.mat');
end
if nargin<3 || isempty(protectRadius), protectRadius=18; end
if nargin<4 || isempty(spacing), spacing=1.5; end
if nargin<5 || isempty(minAngleTarget), minAngleTarget=18; end
assert(protectRadius>0 && spacing>0 && minAngleTarget>0 && minAngleTarget<60);
S=load(sourceFile);
xy=S.pAll; tri=S.tSoil;
originalInterface=S.soilTunnelCoupleTable(:,2:3);
assert(size(xy,2)==2 && size(tri,2)==3 && ...
    all(tri(:)==fix(tri(:))) && all(tri(:)>=1) && ...
    all(tri(:)<=size(xy,1)));
center=[S.xc S.yc];
radius=hypot(xy(:,1)-center(1),xy(:,2)-center(2));
originalBoundary=boundaryEdges(tri);
boundaryNodes=unique(originalBoundary(:));
protect=any(radius(tri)<=protectRadius,2) | ...
    any(ismember(tri,boundaryNodes),2);
protectedTri=tri(protect,:);
assert(~isempty(protectedTri) && ~all(protect));
protectedSoil=unique(protectedTri(:));
soilNodes=unique(tri(:));
fixedSoil=unique([protectedSoil; originalBoundary(:)]);
candidate=setdiff(soilNodes,fixedSoil);
[~,order]=sortrows([radius(candidate) candidate]);
candidate=candidate(order);
retained=fixedSoil(:);
selected=false(numel(candidate),1);
for k=1:numel(candidate)
    node=candidate(k);
    delta=xy(retained,:)-xy(node,:);
    if min(sum(delta.^2,2))>=spacing^2
        retained(end+1,1)=node; %#ok<AGROW>
        selected(k)=true;
    end
end
retained=sort(retained);

% Preserve every edge of each protected triangle as a constraint. The
% protected triangles must appear verbatim in the constrained triangulation.
constraints=unique([originalBoundary;allEdges(protectedTri)],'rows');
soilLocal=zeros(size(xy,1),1);
soilLocal(soilNodes)=1:numel(soilNodes);
TR=triangulation(soilLocal(tri),xy(soilNodes,:));
qualityInsertions=0;
for pass=1:12
    local=zeros(size(xy,1),1); local(retained)=1:numel(retained);
    assert(all(local(constraints(:))>0));
    DT=delaunayTriangulation(xy(retained,:),local(constraints));
    dtTri=retained(DT.ConnectivityList);
    [present,~]=ismember(sort(protectedTri,2),sort(dtTri,2),'rows');
    assert(all(present),'Constrained triangulation changed a protected triangle.');
    centroid=(xy(dtTri(:,1),:)+xy(dtTri(:,2),:)+xy(dtTri(:,3),:))/3;
    inSoil=~isnan(pointLocation(TR,centroid));
    isProtected=ismember(sort(dtTri,2),sort(protectedTri,2),'rows');
    assert(all(inSoil(isProtected)));
    angle=triangleMinAngle(xy,dtTri);
    bad=find(inSoil & ~isProtected & angle<minAngleTarget);
    if isempty(bad), break; end
    remaining=setdiff(soilNodes,retained);
    if isempty(remaining), break; end
    containing=pointLocation(DT,xy(remaining,:));
    add=zeros(numel(bad),1); nadd=0;
    for k=1:numel(bad)
        b=bad(k);
        choices=remaining(containing==b);
        if isempty(choices)
            originalID=pointLocation(TR,centroid(b,:));
            if ~isnan(originalID)
                choices=intersect(tri(originalID,:),remaining);
            end
        end
        if isempty(choices)
            d2=sum((xy(remaining,:)-centroid(b,:)).^2,2);
            [~,j]=min(d2); choices=remaining(j);
        end
        d2=sum((xy(choices,:)-centroid(b,:)).^2,2);
        [~,j]=min(d2);
        nadd=nadd+1; add(nadd)=choices(j);
    end
    add=unique(add(1:nadd));
    add=setdiff(add,retained);
    if isempty(add), break; end
    retained=sort([retained;add]);
    qualityInsertions=qualityInsertions+numel(add);
end
assert(~any(inSoil & ~isProtected & angle<minAngleTarget), ...
    'Quality target could not be met using existing far-field nodes.');
farTri=dtTri(inSoil & ~isProtected,:);
newTriOld=[protectedTri;farTri];
removed=setdiff(soilNodes,retained);
assert(size(unique(sort(newTriOld,2),'rows'),1)==size(newTriOld,1), ...
    'Retrangulation made duplicate triangles.');
assert(isequal(boundaryEdges(newTriOld),originalBoundary), ...
    'Soil-domain boundary edges changed.');
oldArea=triangleArea(xy,tri); newArea=triangleArea(xy,newTriOld);
assert(all(newArea>1e-12));
assert(abs(sum(newArea)-sum(oldArea))<1e-8*sum(oldArea), ...
    'Soil area changed, indicating a gap or overlap.');

% Other mesh blocks contain global node IDs. The interface table also
% records the soil triangle ID (column 4) and its local edge (column 5).
% Protected interface triangles are left in their original ordering.
keepAll=setdiff((1:size(xy,1))',removed);
oldToNew=zeros(size(xy,1),1);
oldToNew(keepAll)=1:numel(keepAll);
S.pAll=xy(keepAll,:);
S.tSoil=oldToNew(newTriOld);
% Coincident soil and lining nodes have different physics. Keep the exact
% original-node IDs; coordinates alone cannot identify an initial field.
S.newNodeToOriginalNode=keepAll;
for field={'qTunnelGlobal','qInvertGlobal','tInvertTriGlobal', ...
        'tunnelInnerRingGlobal','invertBottomRingGlobal', ...
        'invertTopRingGlobal','bInvertTopGlobal'}
    name=field{1};
    if isfield(S,name) && ~isempty(S.(name))
        indices=S.(name);
        assert(all(oldToNew(indices(:))>0), ...
            ['Structural node removed from ' name '.']);
        S.(name)=oldToNew(indices);
    end
end
C=S.soilTunnelCoupleTable;
assert(all(protect(C(:,4))),'A soil-lining interface triangle was not protected.');
oldTriToNew=zeros(size(tri,1),1);
oldTriToNew(find(protect))=1:nnz(protect); %#ok<FNDSB>
C(:,4)=oldTriToNew(C(:,4));
for col=[2 3 6 7]
    assert(all(oldToNew(C(:,col))>0));
    C(:,col)=oldToNew(C(:,col));
end
S.soilTunnelCoupleTable=C;

[negativeOriginal,interfaceWeightOriginal]=negativeHydraulicEdges( ...
    xy,tri,originalInterface);
[negativeCoarse,interfaceWeightCoarse]=negativeHydraulicEdges( ...
    S.pAll,S.tSoil,C(:,2:3));
assert(negativeCoarse.interface<=negativeOriginal.interface, ...
    'Coarsening increased negative interface hydraulic weights.');
assert(max(abs(interfaceWeightOriginal-interfaceWeightCoarse))<1e-10, ...
    'Protected interface hydraulic weights changed.');
report=struct('sourceFile',sourceFile,'outputFile',outputFile, ...
    'protectRadius_m',protectRadius,'minimumFarSpacing_m',spacing, ...
    'minAngleTarget_deg',minAngleTarget, ...
    'originalNodes',size(xy,1),'coarseNodes',size(S.pAll,1), ...
    'originalSoilTriangles',size(tri,1), ...
    'coarseSoilTriangles',size(S.tSoil,1), ...
    'originalSoilArea_m2',sum(oldArea), ...
    'coarseSoilArea_m2',sum(newArea), ...
    'originalBoundaryEdges',size(originalBoundary,1), ...
    'coarseBoundaryEdges',size(boundaryEdges(S.tSoil),1), ...
    'protectedTriangles',nnz(protect), ...
    'farNodesSelected',nnz(selected), ...
    'qualityInsertions',qualityInsertions, ...
    'originalMinimumAngle_deg',min(triangleMinAngle(xy,tri)), ...
    'coarseMinimumAngle_deg',min(triangleMinAngle(xy,newTriOld)), ...
    'negativeHydraulicEdgesOriginal',negativeOriginal, ...
    'negativeHydraulicEdgesCoarse',negativeCoarse);
S.coarseningReport=report;
folder=fileparts(outputFile);
if ~isfolder(folder), mkdir(folder); end
save(outputFile,'-struct','S','-v7');
fprintf('Far-field mesh: nodes %d -> %d; soil triangles %d -> %d; interface negative edges %d -> %d.\n', ...
    report.originalNodes,report.coarseNodes, ...
    report.originalSoilTriangles,report.coarseSoilTriangles, ...
    negativeOriginal.interface,negativeCoarse.interface);
end

function e=allEdges(tri)
e=sort([tri(:,[1 2]);tri(:,[2 3]);tri(:,[3 1])],2);
end

function b=boundaryEdges(tri)
e=allEdges(tri);
[u,~,j]=unique(e,'rows');
count=accumarray(j,1);
b=u(count==1,:);
end

function a=triangleArea(xy,tri)
p=xy(tri(:,1),:); q=xy(tri(:,2),:); r=xy(tri(:,3),:);
cross2=(q(:,1)-p(:,1)).*(r(:,2)-p(:,2)) ...
    -(q(:,2)-p(:,2)).*(r(:,1)-p(:,1));
a=abs(cross2)/2;
end

function angle=triangleMinAngle(xy,tri)
p=xy(tri(:,1),:); q=xy(tri(:,2),:); r=xy(tri(:,3),:);
a=triangleArea(xy,tri);
e12=hypot(p(:,1)-q(:,1),p(:,2)-q(:,2));
e23=hypot(q(:,1)-r(:,1),q(:,2)-r(:,2));
e31=hypot(r(:,1)-p(:,1),r(:,2)-p(:,2));
sinAngle=2*a./[e12.*e31,e12.*e23,e23.*e31];
angle=min(asind(min(1,max(0,sinAngle))),[],2);
end

function [negative,wInterface]=negativeHydraulicEdges(xy,tri,interface)
p=xy(tri(:,1),:); q=xy(tri(:,2),:); r=xy(tri(:,3),:);
a2=(q(:,1)-p(:,1)).*(r(:,2)-p(:,2)) ...
    -(q(:,2)-p(:,2)).*(r(:,1)-p(:,1));
assert(all(abs(a2)>1e-12));
g1=[q(:,2)-r(:,2),r(:,1)-q(:,1)]./a2;
g2=[r(:,2)-p(:,2),p(:,1)-r(:,1)]./a2;
g3=[p(:,2)-q(:,2),q(:,1)-p(:,1)]./a2;
a=abs(a2)/2;
base=[-a.*sum(g1.*g2,2);-a.*sum(g2.*g3,2); ...
    -a.*sum(g3.*g1,2)];
e=allEdges(tri);
[u,~,j]=unique(e,'rows');
w=accumarray(j,base);
negative=struct('all',sum(w< -1e-12),'interface',0);
[found,ix]=ismember(sort(interface,2),u,'rows');
assert(all(found),'Missing soil-lining interface edge.');
wInterface=w(ix);
negative.interface=sum(wInterface< -1e-12);
end
