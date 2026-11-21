function [M,report]=build_segment_mechanical_mesh(S,jointAnglesRad,options)
%BUILD_SEGMENT_MECHANICAL_MESH Eight independent Q4 segment-ring blocks.
% [M,REPORT]=BUILD_SEGMENT_MECHANICAL_MESH(S,ANGLES) creates a mechanical
% mesh only. S, the TH connectivity, and the soil mesh are never changed.
% The eight input joint angles may be unordered; joint IDs retain input order.
% Segment k occupies sortedJointAngle(k) to sortedJointAngle(k+1), wrapping.
%
% Cuts intersect each EXISTING polygonal shell layer, rather than moving its
% vertices onto a new ideal circle. This preserves the original domain and
% gives exact Q4 restriction and affine nodal temperature interpolation.
% M.thermalInterpolation maps original global thermal nodes to mechanical
% ring nodes. Coincident joint nodes have separate mechanical DOFs.
% M.interfaceRemapping maps split outer mechanical edges onto the original
% two-node soil interface; it does not remesh or duplicate the soil.
%
% Only the annular segment ring is included. Invert/slab elements are
% explicitly excluded because they may bridge longitudinal cuts. Requesting
% options.includeInvert=true is rejected instead of silently tying blocks.
% No contact law, bolt, gasket, equilibrium solve or hydraulic feedback is
% introduced by this geometry module. A separate solver must supply them.
if nargin<3, options=struct(); end
assert(isstruct(options)&&isscalar(options), ...
 'SegmentMechanicalMesh:Options','Options must be a scalar struct.');
allowed={'includeInvert','angleTolerance_rad'};
assert(all(ismember(fieldnames(options),allowed)), ...
 'SegmentMechanicalMesh:Options','Unknown mechanical-mesh option.');
includeInvert=false; angleTol=1e-10;
if isfield(options,'includeInvert'), includeInvert=options.includeInvert; end
if isfield(options,'angleTolerance_rad'), angleTol=options.angleTolerance_rad; end
assert(isscalar(includeInvert)&&ismember(includeInvert,[false true]), ...
 'SegmentMechanicalMesh:Options','includeInvert must be a scalar logical.');
assert(~includeInvert,'SegmentMechanicalMesh:InvertUnsupported', ...
 ['This module builds an independent annular ring only. Invert mechanics ' ...
  'require a separate conforming contact/attachment model.']);
assert(isscalar(angleTol)&&isfinite(angleTol)&&angleTol>0&&angleTol<1e-6, ...
 'SegmentMechanicalMesh:Options','Invalid angle tolerance.');
required={'pAll','qTunnelGlobal','xc','yc'};
assert(isstruct(S)&&isscalar(S)&&all(isfield(S,required)), ...
 'SegmentMechanicalMesh:Mesh','A circular-layer Q4 tunnel shell is required.');
xy=double(S.pAll); qOriginal=double(S.qTunnelGlobal); center=double([S.xc S.yc]);
assert(size(xy,2)==2&&all(isfinite(xy(:)))&&numel(center)==2&&all(isfinite(center)), ...
 'SegmentMechanicalMesh:Mesh','Invalid source coordinates or centre.');
validateQuads(qOriginal,size(xy,1));
assert(isnumeric(jointAnglesRad)&&isreal(jointAnglesRad)&& ...
 isvector(jointAnglesRad)&&numel(jointAnglesRad)==8&&all(isfinite(jointAnglesRad(:))), ...
 'SegmentMechanicalMesh:JointAngles','Supply eight finite longitudinal-joint angles.');
inputAngles=mod(double(jointAnglesRad(:)),2*pi);
[joints,sortedToInput]=sort(inputAngles);
jointSpacing=diff([joints;joints(1)+2*pi]);
assert(all(jointSpacing>angleTol), ...
 'SegmentMechanicalMesh:JointAngles','The eight joints must have distinct angular positions.');

[sourceAngles,radii,ring,sourceQuad]=sourceShell(xy,qOriginal,center,angleTol);
nr=numel(radii); nj=8; nOriginal=size(xy,1);
p=zeros(0,2); q=zeros(0,4); nodeSegment=zeros(0,1); elementSegment=zeros(0,1);
nodeRing=zeros(0,1); nodeAngle=zeros(0,1); parentQuad=zeros(0,1);
thermalRows=zeros(0,1); thermalColumns=zeros(0,1); thermalValues=zeros(0,1);
startNodes=zeros(nr,nj); endNodes=zeros(nr,nj);
outerEdges=zeros(0,2); outerParentSector=zeros(0,1); outerSegment=zeros(0,1);
for segment=1:nj
 lo=joints(segment); hi=lo+jointSpacing(segment);
 sourcePhase=lo+mod(sourceAngles-lo,2*pi);
 inside=sort(sourcePhase(sourcePhase>lo+angleTol&sourcePhase<hi-angleTol));
 angles=[lo;inside;hi]; na=numel(angles);
 nodes=reshape(size(p,1)+(1:nr*na),na,nr)';
 for r=1:nr
  for k=1:na
   [point,oldNodes,weights]=ringPoint(angles(k),r,sourceAngles,ring,xy,center,angleTol);
   id=nodes(r,k); p(id,:)=point;
   nodeSegment(id,1)=segment; nodeRing(id,1)=r; nodeAngle(id,1)=mod(angles(k),2*pi);
   thermalRows=[thermalRows;id;id]; %#ok<AGROW>
   thermalColumns=[thermalColumns;oldNodes(:)]; %#ok<AGROW>
   thermalValues=[thermalValues;weights(:)]; %#ok<AGROW>
  end
 end
 for r=1:nr-1
  elements=[nodes(r,1:end-1)',nodes(r,2:end)', ...
   nodes(r+1,2:end)',nodes(r+1,1:end-1)'];
  q=[q;elements]; %#ok<AGROW>
  elementSegment=[elementSegment;segment*ones(na-1,1)]; %#ok<AGROW>
  parents=zeros(na-1,1);
  for k=1:na-1
   sector=sourceSector(.5*(angles(k)+angles(k+1)),sourceAngles,angleTol);
   parents(k)=sourceQuad(r,sector);
  end
  parentQuad=[parentQuad;parents]; %#ok<AGROW>
 end
 startNodes(:,segment)=nodes(:,1); endNodes(:,segment)=nodes(:,end);
 outerEdges=[outerEdges;nodes(1,1:end-1)',nodes(1,2:end)']; %#ok<AGROW>
 outerSegment=[outerSegment;segment*ones(na-1,1)]; %#ok<AGROW>
 sectors=zeros(na-1,1);
 for k=1:na-1
  sectors(k)=sourceSector(.5*(angles(k)+angles(k+1)),sourceAngles,angleTol);
 end
 outerParentSector=[outerParentSector;sectors]; %#ok<AGROW>
end
N=size(p,1);
P=sparse(thermalRows,thermalColumns,thermalValues,N,nOriginal);
assert(all(abs(full(sum(P,2))-1)<1e-12)&&all(nonzeros(P)>=0), ...
 'SegmentMechanicalMesh:Interpolation','Invalid source thermal interpolation.');
assert(norm(P*xy-p,inf)<1e-10, ...
 'SegmentMechanicalMesh:Interpolation','Thermal shape functions do not reproduce coordinates.');

% A joint's minus block precedes it in increasing angle; its plus block
% follows it. Its normal points from minus to plus: e_theta=(-sin,cos).
jointPairs=zeros(nr,2,nj); jointSegments=zeros(nj,2);
interfaceNodes=zeros(nj*(nr-1),4); interfaceID=zeros(nj*(nr-1),1);
interfaceNormal=zeros(nj*(nr-1),2); interfaceTangent=interfaceNormal;
interfaceLength=zeros(nj*(nr-1),1); interfaceSegments=zeros(nj*(nr-1),2);
for sortedID=1:nj
 id=sortedToInput(sortedID); minusSegment=mod(sortedID-2,nj)+1;
 plusSegment=sortedID;
 minus=endNodes(:,minusSegment); plus=startNodes(:,plusSegment);
 jointPairs(:,:,id)=[minus plus]; jointSegments(id,:)=[minusSegment plusSegment];
 rows=(id-1)*(nr-1)+(1:nr-1);
 interfaceNodes(rows,:)=[minus(1:end-1),minus(2:end),plus(1:end-1),plus(2:end)];
 interfaceID(rows)=id;
 interfaceNormal(rows,:)=repmat([-sin(inputAngles(id)),cos(inputAngles(id))],nr-1,1);
 edge=p(minus(2:end),:)-p(minus(1:end-1),:);
 length=hypot(edge(:,1),edge(:,2));
 interfaceLength(rows)=length; interfaceTangent(rows,:)=edge./length;
 interfaceSegments(rows,:)=repmat([minusSegment plusSegment],nr-1,1);
end
minus=reshape(jointPairs(:,1,:),[],1); plus=reshape(jointPairs(:,2,:),[],1);
assert(all(minus~=plus)&&norm(p(minus,:)-p(plus,:),inf)<1e-10, ...
 'SegmentMechanicalMesh:Topology','Joint faces must be coincident with independent nodes.');
assert(all(interfaceLength>0), ...
 'SegmentMechanicalMesh:Topology','A joint interface has zero radial length.');
remapping=outerInterfaceMap(S,p,outerEdges,outerParentSector,outerSegment,ring);
areas=quadAreas(p,q); originalAreas=quadAreas(xy,qOriginal,false);
assert(abs(sum(areas)-sum(originalAreas))<1e-10*max(1,sum(originalAreas)), ...
 'SegmentMechanicalMesh:Geometry','Segment splitting changed the source shell area.');
labels=nodeSegment(q);
assert(isequal(labels,repmat(elementSegment,1,4)), ...
 'SegmentMechanicalMesh:Topology','An element connects two segment blocks.');
components=connectivityComponents(q,N);
assert(numel(unique(components))==nj, ...
 'SegmentMechanicalMesh:Topology','The mechanical ring must have eight disconnected blocks.');
excludedQuad=zeros(0,4); excludedTri=zeros(0,3);
if isfield(S,'qInvertGlobal'), excludedQuad=S.qInvertGlobal; end
if isfield(S,'tInvertTriGlobal'), excludedTri=S.tInvertTriGlobal; end
M=struct('p',p,'q',q,'qTunnelMechanical',q,'nodeSegmentID',nodeSegment, ...
 'elementSegmentID',elementSegment,'nodeRingIndex',nodeRing,'nodeAngleRad',nodeAngle, ...
 'elementParentThermalQ4',parentQuad,'thermalInterpolation',P, ...
 'originalThermalNodeCount',nOriginal,'center',center,'sourceRingRadii_m',radii, ...
 'jointAnglesRad',inputAngles,'sortedJointAnglesRad',joints, ...
 'sortedJointToInputID',sortedToInput,'jointNodePairs',jointPairs, ...
 'jointSegmentPairs',jointSegments, ...
 'jointInterface',struct('nodes',interfaceNodes,'jointID',interfaceID, ...
  'normalMinusToPlus',interfaceNormal,'tangentOuterToInner',interfaceTangent, ...
  'length_m',interfaceLength,'segmentPairs',interfaceSegments), ...
 'interfaceRemapping',remapping,'connectedComponentID',components(:), ...
 'excludedInvertQuadsGlobal',excludedQuad,'excludedInvertTrianglesGlobal',excludedTri, ...
 'scope','Independent annular-ring mechanical geometry only; no invert or soil mechanical elements.', ...
 'constitutiveLawsImplemented',false,'hydraulicFeedbackImplemented',false);
report=struct('segmentCount',nj,'ringLayerCount',nr,'mechanicalNodes',N, ...
 'mechanicalQuads',size(q,1),'originalShellQuads',size(qOriginal,1), ...
 'originalThermalNodes',nOriginal,'disconnectedComponents',numel(unique(components)), ...
 'jointInterfaceElements',size(interfaceNodes,1),'splitOuterEdges',size(outerEdges,1), ...
 'sourceShellArea_m2',sum(originalAreas),'mechanicalShellArea_m2',sum(areas), ...
 'maximumCoordinatePullbackError_m',max(abs(P*xy-p),[],'all'), ...
 'excludedInvertQuadCount',size(excludedQuad,1),'excludedInvertTriangleCount',size(excludedTri,1), ...
 'thermalMeshChanged',false,'soilMeshChanged',false, ...
 'geometryMethod','Ray cuts of existing homothetic polygonal Q4 shell; duplicate mechanical cut nodes', ...
 'limitations','No contact/bolt/gasket constitutive law, mechanical equilibrium, invert mechanics or hydraulic feedback.');
end

function [angles,radii,ring,quadIndex]=sourceShell(xy,quads,center,tol)
ids=unique(quads(:)); relative=xy(ids,:)-center;
r=hypot(relative(:,1),relative(:,2)); a=mod(atan2(relative(:,2),relative(:,1)),2*pi);
a(abs(a-2*pi)<tol)=0;
radii=sort(unique(round(r*1e8)/1e8),'descend'); nr=numel(radii);
assert(nr>=2&&all(radii>0), ...
 'SegmentMechanicalMesh:Mesh','Source shell must have at least two circular vertex layers.');
angles=sort(a(abs(r-radii(1))<1e-7)); na=numel(angles);
assert(na>=8&&all(diff([angles;angles(1)+2*pi])>tol)&& ...
 numel(ids)==nr*na&&size(quads,1)==(nr-1)*na, ...
 'SegmentMechanicalMesh:Mesh','Source shell must be a complete conforming annular Q4 mesh.');
assert(all(diff([angles;angles(1)+2*pi])<pi), ...
 'SegmentMechanicalMesh:Mesh','Shell angular intervals must be smaller than pi.');
ring=zeros(nr,na);
for layer=1:nr
 take=find(abs(r-radii(layer))<1e-7);
 assert(numel(take)==na,'SegmentMechanicalMesh:Mesh','Circular layers differ in angular count.');
 for k=1:na
  match=take(abs(wrapAngle(a(take)-angles(k)))<tol);
  assert(isscalar(match),'SegmentMechanicalMesh:Mesh','Layers do not share angular vertices.');
  ring(layer,k)=ids(match);
 end
end
quadIndex=zeros(nr-1,na); sortedQuads=sort(quads,2);
for layer=1:nr-1
 for k=1:na
  next=mod(k,na)+1;
  expected=sort([ring(layer,k),ring(layer,next),ring(layer+1,next),ring(layer+1,k)]);
  match=find(all(sortedQuads==expected,2));
  assert(isscalar(match),'SegmentMechanicalMesh:Mesh','Missing or duplicated shell Q4 sector.');
  quadIndex(layer,k)=match;
 end
end
end

function [point,ids,weights]=ringPoint(angle,layer,angles,ring,xy,center,tol)
angle=mod(angle,2*pi); match=find(abs(wrapAngle(angles-angle))<tol);
if ~isempty(match)
 assert(isscalar(match)); k=match; next=mod(k,numel(angles))+1;
 weights=[1 0]; ids=ring(layer,[k next]); point=xy(ids(1),:); return
end
k=sourceSector(angle,angles,tol); next=mod(k,numel(angles))+1;
ids=ring(layer,[k next]); A=xy(ids(1),:)-center; edge=xy(ids(2),:)-xy(ids(1),:);
direction=[cos(angle),sin(angle)]; denominator=cross2(direction,edge);
assert(abs(denominator)>1e-14,'SegmentMechanicalMesh:Geometry','A cut is parallel to its ring edge.');
t=cross2(A,direction)/denominator;
assert(t>=-1e-10&&t<=1+1e-10,'SegmentMechanicalMesh:Geometry','A cut missed its ring edge.');
t=max(0,min(1,t)); weights=[1-t t]; point=weights*xy(ids,:);
end

function k=sourceSector(angle,angles,tol)
angle=mod(angle,2*pi); k=find(angles<=angle+tol,1,'last');
if isempty(k), k=numel(angles); end
end

function map=outerInterfaceMap(S,p,edges,parentSector,segment,ring)
ne=size(edges,1); row=zeros(ne,1); soilNodes=zeros(ne,2); concreteNodes=zeros(ne,2);
weights=zeros(ne,2,2); length=zeros(ne,1);
if ~isfield(S,'soilTunnelCoupleTable')||isempty(S.soilTunnelCoupleTable)
 map=struct('available',false,'mechanicalEdgeNodes',edges,'segmentID',segment, ...
  'originalCouplingRow',row,'soilGlobalNodes',soilNodes, ...
  'originalConcreteGlobalNodes',concreteNodes,'soilEndpointWeights',weights,'length_m',length);
 return
end
C=double(S.soilTunnelCoupleTable); xy=double(S.pAll);
assert(size(C,2)>=7,'SegmentMechanicalMesh:Mesh','Invalid soil-lining coupling table.');
concretePair=sort(C(:,6:7),2);
for e=1:ne
 k=parentSector(e); next=mod(k,size(ring,2))+1;
 parent=ring(1,[k next]); match=find(all(concretePair==sort(parent),2));
 assert(isscalar(match),'SegmentMechanicalMesh:Mesh','Outer Q4 edge lacks a unique soil interface.');
 row(e)=match; soilNodes(e,:)=C(match,2:3); concreteNodes(e,:)=C(match,6:7);
 a=xy(soilNodes(e,1),:); edge=xy(soilNodes(e,2),:)-a;
 t=(p(edges(e,:),:)-a)*edge'/sum(edge.^2);
 assert(all(t>=-1e-9&t<=1+1e-9)&& ...
  norm(a+t*edge-p(edges(e,:),:),inf)<1e-9, ...
  'SegmentMechanicalMesh:Mesh','Split mechanical edge is not on its original soil edge.');
 t=max(0,min(1,t)); weights(e,:,:)=reshape([1-t t],1,2,2);
 delta=p(edges(e,2),:)-p(edges(e,1),:); length(e)=hypot(delta(1),delta(2));
end
map=struct('available',true,'mechanicalEdgeNodes',edges,'segmentID',segment, ...
 'originalCouplingRow',row,'soilGlobalNodes',soilNodes, ...
 'originalConcreteGlobalNodes',concreteNodes,'soilEndpointWeights',weights,'length_m',length);
end

function area=quadAreas(xy,quads,requirePositive)
if nargin<3, requirePositive=true; end
x=reshape(xy(quads(:),1),size(quads)); y=reshape(xy(quads(:),2),size(quads));
area=.5*sum(x.*y(:,[2 3 4 1])-y.*x(:,[2 3 4 1]),2);
if requirePositive
 assert(all(area>0),'SegmentMechanicalMesh:Geometry','A mechanical Q4 has nonpositive orientation or area.');
else
 % The legacy TH shell orders its Q4 nodes clockwise and uses |det J|.
 % Reorder new mechanical elements positively, without changing its area.
 assert(all(abs(area)>1e-12),'SegmentMechanicalMesh:Geometry','A source Q4 has zero area.');
 area=abs(area);
end
end

function component=connectivityComponents(quads,N)
a=[quads(:,1);quads(:,2);quads(:,3);quads(:,4)];
b=[quads(:,2);quads(:,3);quads(:,4);quads(:,1)];
component=conncomp(graph(a,b,[],N));
end

function validateQuads(q,N)
assert(isnumeric(q)&&isreal(q)&&size(q,2)==4&&all(isfinite(q(:)))&& ...
 all(q(:)==fix(q(:)))&&all(q(:)>=1&q(:)<=N)&& ...
 all(sum(diff(sort(q,2),1,2)==0,2)==0), ...
 'SegmentMechanicalMesh:Mesh','Invalid Q4 shell connectivity.');
end

function value=cross2(a,b)
value=a(1)*b(2)-a(2)*b(1);
end

function angle=wrapAngle(angle)
angle=mod(angle+pi,2*pi)-pi;
end
