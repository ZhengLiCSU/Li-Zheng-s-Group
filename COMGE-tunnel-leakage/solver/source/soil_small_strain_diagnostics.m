function [d,operators]=soil_small_strain_diagnostics(xy,tri,u)
%SOIL_SMALL_STRAIN_DIAGNOSTICS Complete infinitesimal T3 plane-strain state.
% [D,A]=soil_small_strain_diagnostics(XY,TRI,U) evaluates constant element
% displacement gradients. XY is ns-by-2 [x,y] with y positive DOWN (depth).
% U is ns-by-2 [ux,uy], or [ux1;uy1;ux2;uy2;...] (row vector accepted).
% The caller supplies the displacement INCREMENT relative to the installed
% configuration when checking a THM increment; no reference is subtracted.
%
% D=soil_small_strain_diagnostics(A,U) reuses the sparse operators returned
% by the first call. A.B*U has consecutive [epsxx;epsyy;gammaXY] per element;
% A.gradient*U has [dux/dx;dux/dy;duy/dx;duy/dy] per element. A.V gives div U.
% Orientation reversal of a triangle does not change its physical gradient.
%
% Tensor shear is epsXY=gammaXY/2. Principal strains include the plane-strain
% out-of-plane zero. maximumTensorShearStrain is half the largest 3D
% principal-strain difference; maxAbsTensorShearXY is the shear COMPONENT
% in the supplied coordinates. These are different quantities, e.g. under
% biaxial dilation with constrained out-of-plane strain.
% rotationComponent is (duy/dx-dux/dy)/2; rotationNorm is the Frobenius norm
% of the infinitesimal skew tensor. A pure infinitesimal rigid rotation has
% zero strain but nonzero rotation/gradient. This does NOT validate finite
% rotations, nor change E, boundaries, physical geometry or constitutive laws.

if nargin==2 && isstruct(xy)
    operators=xy; u=tri;
    assert(isscalar(operators)&&all(isfield(operators, ...
        {'B','gradient','V','nodeCount','elementCount'})), ...
        'SoilSmallStrain:Cache','Supply the returned sparse operator cache.');
    ns=operators.nodeCount; ne=operators.elementCount;
    assert(isscalar(ns)&&ns>=3&&ns==fix(ns)&&isscalar(ne)&&ne>=1&&ne==fix(ne) ...
        &&isequal(size(operators.B),[3*ne 2*ns]) ...
        &&isequal(size(operators.gradient),[4*ne 2*ns]) ...
        &&isequal(size(operators.V),[ne 2*ns]), ...
        'SoilSmallStrain:Cache','Invalid sparse operator dimensions.');
else
    assert(nargin==3,'SoilSmallStrain:InputCount', ...
        'Supply XY, TRI, U or the operator cache and U.');
    operators=assembleOperators(xy,tri);
    ns=operators.nodeCount; ne=operators.elementCount;
end
u=displacementVector(u,ns);
strain=reshape(operators.B*u,3,ne)';
gradient=reshape(operators.gradient*u,4,ne)';
exx=strain(:,1); eyy=strain(:,2); gamma=strain(:,3);
center=(exx+eyy)/2;
radius=hypot((exx-eyy)/2,gamma/2);
principal=[center-radius,center+radius,zeros(ne,1)];
elementMaximumPrincipal=max(abs(principal),[],2);
elementMaximumShear=(max(principal,[],2)-min(principal,[],2))/2;
rotation=(gradient(:,3)-gradient(:,2))/2;
gradientNorm=sqrt(sum(gradient.^2,2));
rotationNorm=sqrt(2)*abs(rotation);
volume=exx+eyy;
d=struct('strainComponents',strain,'epsxx',exx,'epsyy',eyy, ...
    'gammaXY',gamma,'tensorShearXY',gamma/2,'volumetricStrain',volume, ...
    'principalStrains',principal, ...
    'elementMaxAbsPrincipalStrain',elementMaximumPrincipal, ...
    'elementMaximumTensorShearStrain',elementMaximumShear, ...
    'displacementGradient',gradient,'gradientNorm',gradientNorm, ...
    'rotationComponent',rotation,'rotationNorm',rotationNorm, ...
    'maxAbsPrincipalStrain',max(elementMaximumPrincipal), ...
    'maximumTensorShearStrain',max(elementMaximumShear), ...
    'maxAbsTensorShearXY',max(abs(gamma/2)), ...
    'maxAbsVolumetricStrain',max(abs(volume)), ...
    'maxGradientNorm',max(gradientNorm),'maxRotationNorm',max(rotationNorm), ...
    'coordinateSystem','x horizontal; y positive-down depth; out-of-plane tunnel-axis strain zero', ...
    'strainColumns',{{'epsxx','epsyy','gammaXY'}}, ...
    'gradientColumns',{{'dux_dx','dux_dy','duy_dx','duy_dy'}}, ...
    'principalColumns',{{'inPlaneMinimum','inPlaneMaximum','outOfPlaneZero'}});
end

function A=assembleOperators(xy,tri)
assert(isnumeric(xy)&&isreal(xy)&&ismatrix(xy)&&size(xy,2)==2 ...
    &&size(xy,1)>=3&&all(isfinite(xy(:))), ...
    'SoilSmallStrain:Mesh','XY must contain finite ns-by-2 coordinates.');
ns=size(xy,1); ne=size(tri,1);
assert(isnumeric(tri)&&isreal(tri)&&ismatrix(tri)&&size(tri,2)==3&&ne>=1 ...
    &&all(isfinite(tri(:)))&&all(tri(:)==fix(tri(:))) ...
    &&all(tri(:)>=1&tri(:)<=ns), ...
    'SoilSmallStrain:Mesh','TRI must contain valid nonempty T3 connectivity.');
xy=double(xy); tri=double(tri);
bi=zeros(18*ne,1); bj=bi; bv=bi;
gi=zeros(24*ne,1); gj=gi; gv=gi;
area=zeros(ne,1); shapeGradients=zeros(ne,2,3);
for e=1:ne
    ids=tri(e,:); x=xy(ids,:);
    cross=(x(2,1)-x(1,1))*(x(3,2)-x(1,2)) ...
        -(x(3,1)-x(1,1))*(x(2,2)-x(1,2));
    assert(isfinite(cross)&&cross~=0,'SoilSmallStrain:Degenerate', ...
        'Each T3 must have a finite nonzero signed area.');
    area(e)=abs(cross)/2;
    g=[x(2,2)-x(3,2),x(3,2)-x(1,2),x(1,2)-x(2,2); ...
        x(3,1)-x(2,1),x(1,1)-x(3,1),x(2,1)-x(1,1)]/cross;
    assert(all(isfinite(g(:))),'SoilSmallStrain:Degenerate', ...
        'Triangle shape gradients must be finite.');
    shapeGradients(e,:,:)=reshape(g,1,2,3);
    gradient=zeros(4,6);
    gradient(1,1:2:6)=g(1,:); gradient(2,1:2:6)=g(2,:);
    gradient(3,2:2:6)=g(1,:); gradient(4,2:2:6)=g(2,:);
    strain=[gradient(1,:);gradient(4,:);gradient(2,:)+gradient(3,:)];
    dof=reshape([2*ids-1;2*ids],[],1);
    rows=repmat((3*e-2:3*e)',1,6); columns=repmat(dof',3,1);
    k=18*(e-1)+(1:18); bi(k)=rows(:); bj(k)=columns(:); bv(k)=strain(:);
    rows=repmat((4*e-3:4*e)',1,6); columns=repmat(dof',4,1);
    k=24*(e-1)+(1:24); gi(k)=rows(:); gj(k)=columns(:); gv(k)=gradient(:);
end
B=sparse(bi,bj,bv,3*ne,2*ns); G=sparse(gi,gj,gv,4*ne,2*ns);
select=sparse(repelem((1:ne)',2),reshape([3*(1:ne)-2;3*(1:ne)-1],[],1),1,ne,3*ne);
A=struct('B',B,'gradient',G,'V',select*B,'nodeCount',ns, ...
    'elementCount',ne,'elementAreas',area,'shapeGradients',shapeGradients, ...
    'strainRowOrder','Consecutive epsxx, epsyy, engineering gammaXY per element', ...
    'gradientRowOrder','Consecutive dux/dx, dux/dy, duy/dx, duy/dy per element');
end

function value=displacementVector(value,ns)
assert(isnumeric(value)&&isreal(value)&&all(isfinite(value(:))) ...
    &&((isvector(value)&&numel(value)==2*ns)||isequal(size(value),[ns 2])), ...
    'SoilSmallStrain:Displacement','U must be finite interleaved 2ns or ns-by-2.');
if isequal(size(value),[ns 2]), value=reshape(value',[],1); end
value=double(value(:));
end
