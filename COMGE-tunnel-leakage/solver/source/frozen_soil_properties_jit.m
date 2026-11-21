function Y=frozen_soil_properties_jit(T,h,v)
%FROZEN_SOIL_PROPERTIES_JIT Numeric, properties-only frozen soil kernel.
%#codegen
% T and h are double column vectors (either may be scalar). v is the fixed
% 22-by-1 scientific-parameter vector from frozen_soil_jit_parameters.
% Y has one row per point and columns
% [sl si sg thetaW thetaI thetaA total elasticWater mass energy k K].
% This reproduces frozen_soil_constitutive(...,false), including equal phase
% densities, Painter-Karra beta=omega=1, strict phase switches, log-domain
% retention, Mualem relative permeability, ice impedance and IAPWS viscosity.
% It changes neither the physics nor SCV integration. There are no tangents,
% tables, fitted substitutes, permeability floors or reduced precision.

assert(isa(T,'double')&&isa(h,'double')&&isa(v,'double'));
assert(isreal(T)&&isreal(h)&&isreal(v));
assert(size(T,2)==1&&size(h,2)==1&&numel(v)==22);
assert(all(isfinite(T(:)))&&all(isfinite(h(:)))&&all(isfinite(v(:))));
assert(all(v([1:5,7:8,12:19])>0)&&all(v([6,9:11])>=0));
assert(v(5)<1&&v(6)<v(5)&&v(8)>1&&(v(20)==0||v(20)==1));
nt=numel(T); nh=numel(h);
assert(nt==nh||nt==1||nh==1);
if nt==1
    n=nh;
elseif nh==1
    n=nt;
else
    n=nt;
end
viscosityEnabled=v(20)==1;
if viscosityEnabled
    assert(all(T(:)>=-20&T(:)<=110)&&v(21)>=-20&&v(21)<=110);
end

rho=v(1); g=v(2); L=v(3); Tm=v(4); phi=v(5); sr=v(6)/phi;
alpha=v(7); vgN=v(8); m=1-1/vgN; Ss=v(9); Ksat=v(10);
iceImpedance=v(11); cw=v(12); ci=v(13); rhoS=v(14); cs=v(15);
kw=v(16); ki=v(17); ks=v(18); ka=v(19); reference=v(22);
solidCapacity=(1-phi)*rhoS*cs;
muRef=1;
if viscosityEnabled, muRef=viscosity(v(21)); end
Y=zeros(n,12);
for i=1:n
    ti=T(min(i,nt)); hi=h(min(i,nh));
    r=retention(hi,alpha,vgN,m,sr);
    freezeHead=L/(g*Tm)*min(ti,0);
    F=retention(freezeHead,alpha,vgN,m,sr);
    sl=r; si=0;
    if F<r
        sl=F;
        si=1-F/r;
    end
    sg=1-sl-si;
    thetaW=phi*sl; thetaI=phi*si; thetaA=phi*sg;
    total=thetaW+thetaI; elastic=Ss*max(hi,0); mass=total+elastic;
    liquidEnthalpy=L+cw*ti; iceEnthalpy=ci*ti;
    energy=solidCapacity*ti+rho*(thetaW*liquidEnthalpy+ ...
        thetaI*iceEnthalpy+elastic*liquidEnthalpy);
    energy=energy+rho*reference*mass;
    k=(1-phi)*ks+thetaW*kw+thetaI*ki+thetaA*ka;
    Se=max(0,min(1,(sl-sr)/(1-sr)));
    b=1-(1-Se^(1/m))^m;
    kr=sqrt(Se)*b^2;
    iceFactor=10^(-iceImpedance*si);
    mobility=1;
    if viscosityEnabled, mobility=muRef/viscosity(ti); end
    K=Ksat*kr*iceFactor*mobility;
    Y(i,:)=[sl si sg thetaW thetaI thetaA total elastic mass energy k K];
end
end

function s=retention(h,alpha,n,m,sr)
s=1;
if h<0
    u=-h; logx=n*(log(alpha)+log(u));
    logOnePlusX=max(logx,0)+log1p(exp(-abs(logx)));
    logq=-m*logOnePlusX;
    s=sr+(1-sr)*exp(logq);
end
end

function mu=viscosity(T)
tk=T+273.15; mu=0;
a=[280.68,511.45,61.131,.45903]; b=[-1.9,-7.7,-19.6,-40];
for j=1:4
    term=1e-6*a(j)*(tk/300)^b(j);
    mu=mu+term;
end
end
