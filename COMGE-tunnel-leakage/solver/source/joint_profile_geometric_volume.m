function [Vgeom,dVdg,details,dVdu]=joint_profile_geometric_volume(B,ur,Cref,fraction)
%JOINT_PROFILE_GEOMETRIC_VOLUME Exact positive-part volume of P1 joint gaps.
% [V,DG,INFO,DU]=joint_profile_geometric_volume(B,UR,CREF,F) takes the cache
% from build_joint_gap_profiles and interleaved ring displacement UR [m].
% Each radial P1 segment is integrated over its ACTUAL nodal-radius length:
%   V_j = f_j * integral max(g_j(r),0) dr      [m3/m = m2].
% The equivalent aperture is f*max(g,0); f is scalar or one value per joint
% in [0,1]. CREF is scalar or nj nonnegative conductances. CREF=0 denotes
% a sealed joint and gives exactly zero volume and derivatives.
%
% DG is nj-by-nr, in the ORIGINAL mechanical radial-node order, with units
% m. DU is sparse nj-by-numel(UR), obtained through B.nodalGapU. It is built
% only if the fourth output is requested. Nodal radii may be ascending or
% descending; positive-part integrals are independent of their orientation.
%
% Both positive endpoints give a trapezoid; both nonpositive endpoints give
% zero. A sign crossing gives L/2*gpositive^2/(gpositive-gnegative). The
% derivative at a zero endpoint next to a nonzero endpoint is continuous;
% at two exactly zero endpoints no unique derivative exists, and the
% inactive (zero) generalized derivative is selected.
%
% Thermal midpoint samples, phase state and wholePathOpen are intentionally
% NOT used. A partially open geometric inventory exists even if every
% thermal sample is closed and the hydraulic through-path is blocked. This
% is a geometry-only ceiling/diagnostic, not a moving-channel storage law.

assert(nargin==4,'JointProfileVolume:InputCount','Supply B, UR, CREF and F.');
assert(isstruct(B)&&isscalar(B)&&all(isfield(B, ...
    {'nr','nj','nodalRadius_m','nodalGap0','nodalGapU'})), ...
    'JointProfileVolume:Cache','Supply the cache returned by build_joint_gap_profiles.');
nr=B.nr; nj=B.nj;
assert(isnumeric(nr)&&isscalar(nr)&&isfinite(nr)&&nr>=2&&nr==fix(nr) ...
    &&isnumeric(nj)&&isscalar(nj)&&isfinite(nj)&&nj>=1&&nj==fix(nj), ...
    'JointProfileVolume:Cache','Invalid radial-node or joint inventory.');
radius=B.nodalRadius_m;
assert(isnumeric(radius)&&isreal(radius)&&isequal(size(radius),[nr nj]) ...
    &&all(isfinite(radius(:)))&&all(radius(:)>0), ...
    'JointProfileVolume:Radius','Nodal radii must be positive finite nr-by-nj values.');
map=B.nodalGapU; gap0=B.nodalGap0;
assert(isnumeric(map)&&isreal(map)&&size(map,1)==nr*nj&&size(map,2)>=1 ...
    &&all(isfinite(nonzeros(map)))&&isnumeric(gap0)&&isreal(gap0) ...
    &&isvector(gap0)&&numel(gap0)==nr*nj&&all(isfinite(gap0(:))), ...
    'JointProfileVolume:Cache','Invalid joint-major nodal gap maps.');
assert(isnumeric(ur)&&isreal(ur)&&isvector(ur)&&numel(ur)==size(map,2) ...
    &&all(isfinite(ur(:))),'JointProfileVolume:Displacement', ...
    'UR must contain all finite ring displacement DOFs.');
conductance=jointValue(Cref,nj,0,Inf,'Conductance');
fraction=jointValue(fraction,nj,0,1,'Fraction');
radius=double(radius); ur=double(ur(:));
nodal=reshape(double(gap0(:))+map*ur,nr,nj)';
assert(all(isfinite(nodal(:))),'JointProfileVolume:Displacement', ...
    'The evaluated nodal gap must remain finite.');
active=conductance>0 & fraction>0;
Vgeom=zeros(nj,1); dVdg=zeros(nj,nr);
segmentLength=zeros(nj,nr-1); segmentVolume=segmentLength;
radialOrder=zeros(nj,nr); actualLength=zeros(nj,1);
for joint=1:nj
    [r,order]=sort(radius(:,joint));
    length=diff(r);
    assert(all(length>0),'JointProfileVolume:Radius', ...
        'Every joint must have distinct mechanical radial nodes.');
    radialOrder(joint,:)=order'; segmentLength(joint,:)=length';
    actualLength(joint)=r(end)-r(1);
    if ~active(joint), continue; end
    for k=1:nr-1
        a=order(k); b=order(k+1); ga=nodal(joint,a); gb=nodal(joint,b);
        halfLength=length(k)/2;
        integral=0; da=0; db=0;
        if ga>=0 && gb>=0 && (ga>0 || gb>0)
            integral=length(k)*(.5*ga+.5*gb);
            da=halfLength; db=halfLength;
        elseif ga>0 && gb<0
            ratio=positiveFraction(ga,gb);
            integral=halfLength*ga*ratio;
            da=halfLength*ratio*(2-ratio); db=halfLength*ratio^2;
        elseif gb>0 && ga<0
            ratio=positiveFraction(gb,ga);
            integral=halfLength*gb*ratio;
            db=halfLength*ratio*(2-ratio); da=halfLength*ratio^2;
        end
        segmentVolume(joint,k)=fraction(joint)*integral;
        Vgeom(joint)=Vgeom(joint)+segmentVolume(joint,k);
        dVdg(joint,a)=dVdg(joint,a)+fraction(joint)*da;
        dVdg(joint,b)=dVdg(joint,b)+fraction(joint)*db;
    end
end
assert(all(isfinite(Vgeom))&&all(isfinite(dVdg(:))), ...
    'JointProfileVolume:Overflow','Geometric volume or its derivatives overflow.');
aperture=fraction.*max(nodal,0).*double(conductance>0);
details=struct('nodalGap_m',nodal,'nodalAperture_m',aperture, ...
    'nodalRadius_m',radius','ascendingToOriginalRadialOrder',radialOrder, ...
    'segmentLength_m',segmentLength,'segmentVolume_m2',segmentVolume, ...
    'actualPathLength_m',actualLength,'geometricVolume_m2',Vgeom, ...
    'dVolume_dNodalGap_m',dVdg,'sealedMask',conductance==0, ...
    'openingTransmissionFraction',fraction, ...
    'wholePathOpen',all(nodal>0,2), ...
    'integrationMethod','Exact positive-part integral of each physical radial P1 segment; independent of thermal samples', ...
    'nodalOrder','nj-by-nr, original mechanical radial-node order; segment columns follow ascending radius', ...
    'zeroGapDerivative','Inactive generalized derivative when both segment endpoints are exactly zero');
dVdu=[];
if nargout>=4
    rows=repelem((1:nj)',nr); columns=(1:nj*nr)';
    nodalDerivative=sparse(rows,columns,reshape(dVdg',[],1),nj,nj*nr);
    dVdu=nodalDerivative*sparse(map);
end
end

function ratio=positiveFraction(positive,negative)
% Scaling avoids overflow in positive-negative while preserving its ratio.
scale=max(positive,-negative);
ratio=(positive/scale)/((positive/scale)+(-negative/scale));
end

function value=jointValue(value,nj,minimum,maximum,name)
assert(isnumeric(value)&&isreal(value)&&isvector(value)&&all(isfinite(value(:))) ...
    &&all(value(:)>=minimum&value(:)<=maximum)&&(isscalar(value)||numel(value)==nj), ...
    ['JointProfileVolume:' name],'Invalid scalar or per-joint %s.',name);
value=double(value(:));
if isscalar(value), value=repmat(value,nj,1); end
end
