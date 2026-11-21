function G = build_soil_cvfe_flux(xy,tri,soilNodes)
%BUILD_SOIL_CVFE_FLUX Barycentric median-dual faces for triangular P1 fields.
% G=build_soil_cvfe_flux(xy,tri,soilNodes) preserves the supplied local soil
% numbering. xy and tri use global numbering, including any separate lining
% nodes. Each triangle contributes three distinct internal SCV subfaces,
% from its edge midpoint to its centroid. These are NOT combined across
% triangles: a face-normal Darcy flux must remain distinct during upwinding.
%
% G.D0 maps local total hydraulic potential H into physical subface flows
% for unit conductivity. G.inc deposits outward flux at the originating
% node and its negative at the receiving node. Thus G.inc*D(K)*H equals the
% P1 Galerkin diffusion residual for elementwise scalar K, including on
% obtuse triangles. Flux has units of area/time per unit axial tunnel length.
% Coordinates follow the existing tunnel convention: y is positive down,
% so pressure head h gives total hydraulic potential H=h-y.

validateattributes(xy,{'numeric'},{'2d','real','finite','nonempty'}, ...
    mfilename,'xy',1);
validateattributes(tri,{'numeric'},{'2d','real','finite','integer','positive', ...
    'nonempty'},mfilename,'tri',2);
if size(xy,2)~=2 || size(tri,2)~=3 || any(tri(:)>size(xy,1))
    error('build_soil_cvfe_flux:Geometry', ...
        'xy must be N-by-2 and tri must contain three valid global node IDs.');
end
N=size(xy,1); ne=size(tri,1);
if nargin<3 || isempty(soilNodes)
    soilNodes=unique(tri(:));
end
validateattributes(soilNodes,{'numeric'},{'vector','real','finite','integer', ...
    'positive','nonempty'},mfilename,'soilNodes',3);
soilNodes=soilNodes(:); ns=numel(soilNodes);
if any(soilNodes>N) || numel(unique(soilNodes))~=ns
    error('build_soil_cvfe_flux:SoilNodes', ...
        'soilNodes must be unique valid global node IDs.');
end
globalToLocal=zeros(N,1); globalToLocal(soilNodes)=(1:ns)';
localTri=reshape(globalToLocal(tri(:)),ne,3);
if any(localTri(:)==0)
    error('build_soil_cvfe_flux:SoilNodes', ...
        'Every triangle node must occur in soilNodes.');
end
x=reshape(xy(tri(:),1),ne,3); y=reshape(xy(tri(:),2),ne,3);
twiceSignedArea=(x(:,2)-x(:,1)).*(y(:,3)-y(:,1)) ...
    -(x(:,3)-x(:,1)).*(y(:,2)-y(:,1));
edgeScale=max([(x(:,2)-x(:,1)).^2+(y(:,2)-y(:,1)).^2, ...
    (x(:,3)-x(:,2)).^2+(y(:,3)-y(:,2)).^2, ...
    (x(:,1)-x(:,3)).^2+(y(:,1)-y(:,3)).^2],[],2);
if any(abs(twiceSignedArea)<=64*eps.*edgeScale)
    error('build_soil_cvfe_flux:DegenerateTriangle', ...
        'The mesh contains a degenerate or numerically singular triangle.');
end
area=abs(twiceSignedArea)/2;
dNx=[y(:,2)-y(:,3),y(:,3)-y(:,1),y(:,1)-y(:,2)]./twiceSignedArea;
dNy=[x(:,3)-x(:,2),x(:,1)-x(:,3),x(:,2)-x(:,1)]./twiceSignedArea;
centroid=[mean(x,2),mean(y,2)];
pairs=[1 2;2 3;3 1]; nf=3*ne;
from=zeros(nf,1); to=from; normalLength=zeros(nf,2);
edgeMidpoint=normalLength; sampleXY=normalLength;
for pair=1:3
    faces=(pair-1)*ne+(1:ne)'; i=pairs(pair,1); j=pairs(pair,2);
    midpoint=[(x(:,i)+x(:,j))/2,(y(:,i)+y(:,j))/2];
    subfaceVector=centroid-midpoint;
    normal=[subfaceVector(:,2),-subfaceVector(:,1)];
    ijVector=[x(:,j)-x(:,i),y(:,j)-y(:,i)];
    flip=sum(normal.*ijVector,2)<0;
    normal(flip,:)=-normal(flip,:);
    from(faces)=localTri(:,i); to(faces)=localTri(:,j);
    normalLength(faces,:)=normal;
    edgeMidpoint(faces,:)=midpoint;
    sampleXY(faces,:)=(midpoint+centroid)/2;
end
elementIndex=repmat((1:ne)',3,1); faceID=(1:nf)';
faceTri=repmat(localTri,3,1);
coeff=-normalLength(:,1).*repmat(dNx,3,1) ...
    -normalLength(:,2).*repmat(dNy,3,1);
D0=sparse(repmat(faceID,3,1),faceTri(:),coeff(:),nf,ns);
inc=sparse([from;to],[faceID;faceID],[ones(nf,1);-ones(nf,1)],ns,nf);
meanK=sparse(repmat((1:ne)',3,1),localTri(:),ones(3*ne,1)/3,ne,ns);
elementSelect=sparse(faceID,elementIndex,ones(nf,1),nf,ne);
soilSelect=sparse((1:ns)',soilNodes,ones(ns,1),ns,N);
volume=accumarray(localTri(:),repmat(area/3,3,1),[ns 1]);
G=struct('N',N,'ns',ns,'ne',ne,'nf',nf,'soilNodes',soilNodes, ...
    'triGlobal',tri,'localTri',localTri,'xy',xy(soilNodes,:), ...
    'z',xy(soilNodes,2),'area',area,'volume',volume, ...
    'gradX',dNx,'gradY',dNy,'centroid',centroid, ...
    'elementIndex',elementIndex,'from',from,'to',to, ...
    'fromGlobal',soilNodes(from),'toGlobal',soilNodes(to), ...
    'normalLength',normalLength,'edgeMidpoint',edgeMidpoint, ...
    'sampleXY',sampleXY,'inc',inc,'soilSelect',soilSelect, ...
    'meanK',meanK,'elementSelect',elementSelect,'D0',D0);
end
