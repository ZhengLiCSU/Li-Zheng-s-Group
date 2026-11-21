function a=frozen_soil_constitutive(T,h,p,wantDerivatives)
%FROZEN_SOIL_CONSTITUTIVE Equal-density TH properties at arbitrary points.
% T [deg C] and liquid pressure head h [m] have matching sizes; either may
% be scalar. Outputs are densities per bulk volume, WITHOUT mesh weights.
% mass [m3 water-equivalent/m3 bulk], energy [J/m3], k [W/(m K)], K [m/s].
% Painter & Karra (2014), Eqs 15 and 18-19, beta=omega=1:
%   r=R(h), F=R(min(T,0)*L/(g*Tm)), sl=min(r,F), si=1-sl/r.
% The thermal retention argument must NOT also contain h: this explicit
% special case already accounts for pre-freezing retention through min(r,F).
% https://doi.org/10.2136/vzj2013.04.0071
% Water and ice have the SAME density p.rho here. This is the v14 TH closure;
% no mechanics, pressure-melting correction or frost damage is introduced.
% Analytic derivatives are pointwise partials at fixed other primary field.
% At switches choose the nonfrozen min branch (F>=r), the flat retention
% branch (h>=0), and inactive elastic/thermal branches (h<=0, T>=0).
% Away from switches these are classical derivatives. Mualem's
% derivative need not be bounded on approach to full saturation; no claim
% of global smoothness or artificial smoothing of the retention law is made.
% Optional fourth logical input defaults to true. With false, the same
% properties are evaluated without tangents; derivative fields are empty.

if nargin<4, wantDerivatives=true; end
assert(islogical(wantDerivatives)&&isscalar(wantDerivatives), ...
    'frozen_soil_constitutive:DerivativeSwitch','wantDerivatives must be a scalar logical.');
assert(isnumeric(T)&&isreal(T)&&all(isfinite(T(:)))&& ...
    isnumeric(h)&&isreal(h)&&all(isfinite(h(:))), ...
    'frozen_soil_constitutive:InvalidState','T and h must be finite real arrays.');
if isscalar(T)&&~isscalar(h), T=repmat(T,size(h)); end
if isscalar(h)&&~isscalar(T), h=repmat(h,size(T)); end
assert(isequal(size(T),size(h)), ...
    'frozen_soil_constitutive:SizeMismatch','T and h must have matching sizes or be scalar.');
T=double(T); h=double(h); validateParameters(p);
phi=p.porosity; sr=p.thetaR/phi; m=1-1/p.vgN;
freezeHead=p.L/(p.g*p.Tm)*min(T,0);
if wantDerivatives
    [r,rH]=retentionWithDerivative(h,p.vgAlpha,p.vgN,m,sr);
    [F,Fhead]=retentionWithDerivative(freezeHead,p.vgAlpha,p.vgN,m,sr);
else
    r=retentionWithDerivative(h,p.vgAlpha,p.vgN,m,sr);
    F=retentionWithDerivative(freezeHead,p.vgAlpha,p.vgN,m,sr);
end
frozen=F<r;
sl=r; sl(frozen)=F(frozen);
si=zeros(size(T)); si(frozen)=1-F(frozen)./r(frozen);
sg=1-sl-si;

a.sl=sl; a.si=si; a.sg=sg;
a.thetaW=phi*sl; a.thetaI=phi*si; a.thetaA=phi*sg;
a.total=a.thetaW+a.thetaI;
a.elasticWater=p.Ss*max(h,0);
a.mass=a.total+a.elasticWater;
solidCapacity=(1-phi)*p.rhoS*p.cs;
liquidEnthalpy=p.L+p.cw*T; iceEnthalpy=p.ci*T;
a.energy=solidCapacity*T+p.rho*(a.thetaW.*liquidEnthalpy+ ...
    a.thetaI.*iceEnthalpy+a.elasticWater.*liquidEnthalpy);
reference=0;
if isfield(p,'enthalpyReference_Jkg'), reference=p.enthalpyReference_Jkg; end
a.energy=a.energy+p.rho*reference*a.mass;
a.k=(1-phi)*p.ks+a.thetaW*p.kw+a.thetaI*p.ki+a.thetaA*p.ka;

Se=max(0,min(1,(sl-sr)/(1-sr)));
b=1-(1-Se.^(1/m)).^m;
kr=sqrt(Se).*b.^2;
iceFactor=10.^(-p.iceImpedance*si);
mobility=ones(size(T));
if isfield(p,'soilViscosity')&&p.soilViscosity.enabled
    muRef=water_viscosity_iapws(p.soilViscosity.referenceTemperature_C);
    if wantDerivatives
        [mu,dmu]=water_viscosity_iapws(T);
    else
        mu=water_viscosity_iapws(T);
    end
    mobility=muRef./mu;
end
a.K=p.Ksat*kr.*iceFactor.*mobility;
if ~wantDerivatives
    derivatives={'slT','slH','siT','siH','dMassT','dMassH', ...
        'dEnergyT','dEnergyH','dkT','dkH','dKT','dKH'};
    for j=1:numel(derivatives), a.(derivatives{j})=[]; end
    return
end
FT=Fhead*(p.L/(p.g*p.Tm)); FT(T>=0)=0;
slT=zeros(size(T)); slT(frozen)=FT(frozen);
slH=rH; slH(frozen)=0;
siT=zeros(size(T)); siT(frozen)=-slT(frozen)./r(frozen);
siH=zeros(size(T)); siH(frozen)=sl(frozen).*rH(frozen)./r(frozen).^2;
a.slT=slT; a.slH=slH; a.siT=siT; a.siH=siH;
elasticH=p.Ss*double(h>0); % zero is the flat branch at h=0
a.dMassT=phi*(slT+siT);
a.dMassH=phi*(slH+siH)+elasticH;
a.dEnergyT=solidCapacity+p.rho*( ...
    phi*(slT.*liquidEnthalpy+siT.*iceEnthalpy+p.cw*sl+p.ci*si)+ ...
    p.cw*a.elasticWater);
a.dEnergyH=p.rho*(phi*(slH.*liquidEnthalpy+siH.*iceEnthalpy)+ ...
    elasticH.*liquidEnthalpy);
a.dEnergyT=a.dEnergyT+p.rho*reference*a.dMassT;
a.dEnergyH=a.dEnergyH+p.rho*reference*a.dMassH;
a.dkT=phi*((p.kw-p.ka)*slT+(p.ki-p.ka)*siT);
a.dkH=phi*((p.kw-p.ka)*slH+(p.ki-p.ka)*siH);
dkr=zeros(size(T)); inside=Se>0&Se<1;
s=Se(inside); bs=b(inside);
db=(1-s.^(1/m)).^(m-1).*s.^(1/m-1);
dkr(inside)=.5*bs.^2./sqrt(s)+2*sqrt(s).*bs.*db;
SeT=slT/(1-sr); SeH=slH/(1-sr);
SeT(~inside)=0; SeH(~inside)=0;
mobilityT=zeros(size(T));
if isfield(p,'soilViscosity')&&p.soilViscosity.enabled
    mobilityT=-mobility.*dmu./mu;
end
a.dKT=p.Ksat*iceFactor.*(mobility.*(dkr.*SeT- ...
    log(10)*p.iceImpedance*kr.*siT)+kr.*mobilityT);
a.dKH=p.Ksat*iceFactor.*mobility.*(dkr.*SeH- ...
    log(10)*p.iceImpedance*kr.*siH);
end

function [s,ds]=retentionWithDerivative(h,alpha,n,m,sr)
% Log-domain evaluation also supports finite heads whose h^n overflows.
s=ones(size(h)); negative=h<0;
u=-h(negative); logx=n*(log(alpha)+log(u));
logOnePlusX=max(logx,0)+log1p(exp(-abs(logx)));
logq=-m*logOnePlusX;
s(negative)=sr+(1-sr)*exp(logq);
if nargout>1
    ds=zeros(size(h));
    ds(negative)=exp(log((1-sr)*m*n)+logq+ ...
        logx-logOnePlusX-log(u));
end
end

function validateParameters(p)
assert(isstruct(p)&&isscalar(p),'frozen_soil_constitutive:InvalidParameters', ...
    'p must be a scalar parameter structure.');
positive={'rho','g','L','Tm','porosity','vgAlpha','vgN','cw','ci', ...
    'rhoS','cs','kw','ki','ks','ka'};
nonnegative={'thetaR','Ss','Ksat','iceImpedance'};
for names={positive,nonnegative}
    set=names{1};
    for j=1:numel(set)
        name=set{j};
        assert(isfield(p,name)&&isnumeric(p.(name))&&isreal(p.(name))&& ...
            isscalar(p.(name))&&isfinite(p.(name)), ...
            'frozen_soil_constitutive:InvalidParameters', ...
            'Parameter %s must be a finite real scalar.',name);
        if isequal(set,positive), valid=p.(name)>0; else, valid=p.(name)>=0; end
        assert(valid,'frozen_soil_constitutive:InvalidParameters', ...
            'Parameter %s has an invalid sign.',name);
    end
end
assert(p.porosity<1&&p.thetaR<p.porosity&&p.vgN>1, ...
    'frozen_soil_constitutive:InvalidParameters', ...
    'Require 0<porosity<1, 0<=thetaR<porosity, and vgN>1.');
if isfield(p,'enthalpyReference_Jkg')
    assert(isnumeric(p.enthalpyReference_Jkg)&&isreal(p.enthalpyReference_Jkg)&& ...
        isscalar(p.enthalpyReference_Jkg)&&isfinite(p.enthalpyReference_Jkg), ...
        'frozen_soil_constitutive:InvalidParameters', ...
        'enthalpyReference_Jkg must be a finite real scalar.');
end
if isfield(p,'soilViscosity')
    v=p.soilViscosity;
    assert(isstruct(v)&&isscalar(v)&&isfield(v,'enabled')&& ...
        isscalar(v.enabled)&&any(v.enabled==[0 1]), ...
        'frozen_soil_constitutive:InvalidParameters', ...
        'soilViscosity.enabled must be a scalar logical switch.');
    if v.enabled
        assert(isfield(v,'referenceTemperature_C')&& ...
            isnumeric(v.referenceTemperature_C)&&isreal(v.referenceTemperature_C)&& ...
            isscalar(v.referenceTemperature_C)&&isfinite(v.referenceTemperature_C), ...
            'frozen_soil_constitutive:InvalidParameters', ...
            'The enabled viscosity correction requires referenceTemperature_C.');
    end
end
end
