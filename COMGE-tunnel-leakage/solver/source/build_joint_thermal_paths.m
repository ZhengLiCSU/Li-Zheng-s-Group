function paths=build_joint_thermal_paths(S,anglesRad,nq)
%BUILD_JOINT_THERMAL_PATHS Sample concrete temperature through each joint.
% paths.interpolation*T contains nq radial midpoint temperatures per joint,
% with all samples for joint 1 followed by all samples for joint 2, etc.
% Coordinates use the mesh convention; no assumption about the sign of Z
% is made. Angles are measured by ATAN2 about [S.xc,S.yc].
%
% The radial limits follow the actual polygonal Q4 lining, rather than
% the nominal circular radii. qTunnelGlobal identifies the segment shell;
% an invert slab inside that shell is not part of the joint thickness.
% If this shell set is absent, the outermost continuous radial interval
% in all supplied concrete elements is used. Interpolation uses Q4 inverse
% isoparametric coordinates or T3 barycentric coordinates. Soil nodes and
% nearest-node substitutes are never used.
assert(isstruct(S)&&isscalar(S),'Mesh must be a scalar struct.');
assert(isfield(S,'pAll')&&isnumeric(S.pAll)&&isreal(S.pAll) ...
    &&size(S.pAll,1)>0&&size(S.pAll,2)==2&&all(isfinite(S.pAll(:))), ...
    'Mesh must supply finite pAll(:,1:2).');
assert(isfield(S,'xc')&&isfield(S,'yc')&&isreal(S.xc)&&isreal(S.yc) ...
    &&isscalar(S.xc)&&isscalar(S.yc) ...
    &&isfinite(S.xc)&&isfinite(S.yc),'Mesh must supply xc and yc.');
assert(isvector(anglesRad)||isempty(anglesRad),'Joint angles must be a vector.');
assert(isnumeric(anglesRad)&&isreal(anglesRad)&&all(isfinite(anglesRad(:))), ...
    'Joint angles must be real and finite.');
assert(isscalar(nq)&&isfinite(nq)&&nq>=1&&nq==fix(nq), ...
    'nq must be a positive integer.');
xy=double(S.pAll); center=double([S.xc S.yc]); N=size(xy,1);
anglesRad=double(anglesRad(:)); nj=numel(anglesRad); np=nj*nq;
quads=zeros(0,4); triangles=zeros(0,3); shellQuads=zeros(0,4);
for name={'qTunnelGlobal','qInvertGlobal'}
    if isfield(S,name{1})&&~isempty(S.(name{1}))
        ids=double(S.(name{1}));
        validateElements(ids,4,N,name{1});
        quads=[quads;ids]; %#ok<AGROW>
        if strcmp(name{1},'qTunnelGlobal'), shellQuads=ids; end
    end
end
if isfield(S,'tInvertTriGlobal')&&~isempty(S.tInvertTriGlobal)
    triangles=double(S.tInvertTriGlobal);
    validateElements(triangles,3,N,'tInvertTriGlobal');
end
assert(~isempty(quads)||~isempty(triangles)||nj==0, ...
    'No concrete elements are available for thermal path interpolation.');
scale=max(1,max(abs(xy(:)))); coordTol=1e-10*scale;
if ~isempty(shellQuads)
    pathQuads=shellQuads; pathTriangles=zeros(0,3);
else
    pathQuads=quads; pathTriangles=triangles;
end
edges=concreteBoundaryEdges(pathQuads,pathTriangles);
allBounds=elementBounds(xy,quads,triangles);
pathBounds=elementBounds(xy,pathQuads,pathTriangles);
outer=zeros(nj,1); inner=zeros(nj,1); sampleXY=zeros(np,2);
row=zeros(4*np,1); column=row; value=row; used=0;
for joint=1:nj
    direction=[cos(anglesRad(joint)),sin(anglesRad(joint))];
    radii=rayBoundaryIntersections(center,direction,xy,edges,coordTol);
    % Interior edges from a nonconforming concrete mesh can remain in the
    % edge inventory. Membership tests merge their adjacent solid intervals.
    intervals=zeros(0,2);
    for k=1:numel(radii)-1
        if radii(k+1)-radii(k)<=coordTol, continue; end
        point=center+.5*(radii(k)+radii(k+1))*direction;
        [ids,~]=containingConcrete(point,xy,pathQuads,pathTriangles, ...
            pathBounds,coordTol);
        if isempty(ids), continue; end
        if ~isempty(intervals)&&abs(intervals(end,2)-radii(k))<=coordTol
            intervals(end,2)=radii(k+1);
        else
            intervals(end+1,:)=[radii(k),radii(k+1)]; %#ok<AGROW>
        end
    end
    assert(~isempty(intervals), ...
        'Joint %d has no concrete radial path at angle %.9g rad.',joint,anglesRad(joint));
    [~,last]=max(intervals(:,2)); inner(joint)=intervals(last,1);
    outer(joint)=intervals(last,2);
    assert(inner(joint)>coordTol&&outer(joint)-inner(joint)>coordTol, ...
        'Joint %d does not cross an annular concrete segment shell.',joint);
    for k=1:nq
        index=(joint-1)*nq+k;
        % Sample from outside toward the tunnel, in flow-path order.
        radius=outer(joint)-(k-.5)/nq*(outer(joint)-inner(joint));
        point=center+radius*direction; sampleXY(index,:)=point;
        [ids,weights]=containingConcrete(point,xy,quads,triangles, ...
            allBounds,coordTol);
        assert(~isempty(ids), ...
            'Joint %d radial sample %d is outside the concrete mesh.',joint,k);
        assert(all(weights>=-1e-12)&&abs(sum(weights)-1)<1e-10, ...
            'Invalid concrete shape weights at joint %d, sample %d.',joint,k);
        weights=max(weights,0); weights=weights/sum(weights);
        assert(norm(weights*xy(ids,:)-point,inf)<=5*coordTol, ...
            'Concrete interpolation does not reproduce its sample coordinates.');
        slots=used+(1:numel(ids)); row(slots)=index; column(slots)=ids;
        value(slots)=weights; used=used+numel(ids);
    end
end
P=sparse(row(1:used),column(1:used),value(1:used),np,N);
assert(all(abs(full(sum(P,2))-1)<1e-10)&&all(nonzeros(P)>=0), ...
    'Joint thermal interpolation must be a convex combination of concrete DOFs.');
paths=struct('interpolation',P,'sampleXY',sampleXY,'nj',nj,'nq',nq, ...
    'radiusOuter_m',outer,'radiusInner_m',inner);
end

function validateElements(ids,nnode,N,name)
assert(size(ids,2)==nnode&&isreal(ids)&&all(isfinite(ids(:))) ...
    &&all(ids(:)==fix(ids(:)))&&all(ids(:)>=1&ids(:)<=N), ...
    'Invalid %s connectivity.',name);
assert(all(sum(diff(sort(ids,2),1,2)==0,2)==0), ...
    '%s contains an element with repeated nodes.',name);
end

function edges=concreteBoundaryEdges(quads,triangles)
edges=[quads(:,[1 2]);quads(:,[2 3]);quads(:,[3 4]);quads(:,[4 1]); ...
    triangles(:,[1 2]);triangles(:,[2 3]);triangles(:,[3 1])];
if isempty(edges), return; end
[uniqueEdges,~,group]=unique(sort(edges,2),'rows');
count=accumarray(group,1);
assert(all(count<=2),'Concrete connectivity contains a nonmanifold edge.');
edges=uniqueEdges(count==1,:);
end

function bounds=elementBounds(xy,quads,triangles)
bounds.quad=zeros(size(quads,1),4); bounds.triangle=zeros(size(triangles,1),4);
for k=1:size(quads,1)
    x=xy(quads(k,:),:); bounds.quad(k,:)=[min(x,[],1),max(x,[],1)];
end
for k=1:size(triangles,1)
    x=xy(triangles(k,:),:); bounds.triangle(k,:)=[min(x,[],1),max(x,[],1)];
end
end

function radii=rayBoundaryIntersections(center,direction,xy,edges,tol)
radii=0;
for k=1:size(edges,1)
    a=xy(edges(k,1),:)-center; b=xy(edges(k,2),:)-center; v=b-a;
    determinant=cross2(direction,v);
    if abs(determinant)>tol
        radius=cross2(a,v)/determinant;
        edgeFraction=cross2(a,direction)/determinant;
        edgeTol=tol/max(norm(v),tol);
        if radius>=-tol&&edgeFraction>=-edgeTol&&edgeFraction<=1+edgeTol
            radii(end+1,1)=max(radius,0); %#ok<AGROW>
        end
    elseif abs(cross2(a,direction))<=tol
        % A ray along a polygon edge still has its exact endpoint crossings.
        projections=[dot(a,direction);dot(b,direction)];
        radii=[radii;max(projections(projections>=-tol),0)]; %#ok<AGROW>
    end
end
radii=sort(radii(:));
if isempty(radii), radii=0; return; end
keep=[true;diff(radii)>tol]; radii=radii(keep);
end

function [ids,weights]=containingConcrete(point,xy,quads,triangles,bounds,tol)
ids=[]; weights=[];
q=find(inBounds(point,bounds.quad,tol));
for k=q(:)'
    x=xy(quads(k,:),:); uv=[0;0]; ok=false;
    for iter=1:30
        [shape,gradient]=quadShape(uv(1),uv(2));
        errorVector=(point-shape*x)';
        if norm(errorVector,inf)<=tol, ok=true; break; end
        J=(gradient*x)';
        if rcond(J)<1e-13, break; end
        uv=uv+J\errorVector;
        if any(~isfinite(uv))||max(abs(uv))>1e4, break; end
    end
    if ok&&all(abs(uv)<=1+1e-10)
        uv=max(-1,min(1,uv)); [shape,~]=quadShape(uv(1),uv(2));
        if norm(shape*x-point,inf)<=5*tol
            ids=quads(k,:); weights=shape; return;
        end
    end
end
t=find(inBounds(point,bounds.triangle,tol));
for k=t(:)'
    x=xy(triangles(k,:),:); A=[(x(2,:)-x(1,:))',(x(3,:)-x(1,:))'];
    if rcond(A)<1e-13, continue; end
    uv=A\(point-x(1,:))'; shape=[1-sum(uv),uv'];
    if all(shape>=-1e-10)&&all(shape<=1+1e-10)
        shape=max(shape,0); shape=shape/sum(shape);
        if norm(shape*x-point,inf)<=5*tol
            ids=triangles(k,:); weights=shape; return;
        end
    end
end
end

function yes=inBounds(point,bounds,tol)
yes=point(1)>=bounds(:,1)-tol&point(2)>=bounds(:,2)-tol ...
    &point(1)<=bounds(:,3)+tol&point(2)<=bounds(:,4)+tol;
end

function [N,d]=quadShape(xi,eta)
N=.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
d=.25*[-(1-eta),1-eta,1+eta,-(1+eta); ...
    -(1-xi),-(1+xi),1+xi,1-xi];
end

function value=cross2(a,b)
value=a(1)*b(2)-a(2)*b(1);
end
