function [s,d] = cohesive_threephase_state(z,model,p,history,loadingAtHistory)
%THREEPHASE_INTERFACE_STATE Condensed Q4 mechanics + conservative phase EOS.
% z=[log(p_abs/p_atm); ice fraction of CONDENSED mass; vapour VOLUME fraction;
% joint T; solid T]. Vapour mass is retained, including sublimated ice.
if nargin<5, loadingAtHistory=[]; end
n=numel(model.weight); ns=numel(model.capacity); nz=4*n+ns;
assert(numel(z)==nz,'ThreePhase:Shape','Incorrect state length.');
il=1:n; iff=n+(1:n); iv=2*n+(1:n); it=3*n+(1:n); its=4*n+(1:ns);
pa=p.atmosphericPressure_Pa*exp(z(il)); pg=pa-p.atmosphericPressure_Pa;
f=z(iff); theta=z(iv); T=z(it); Ts=z(its);
assert(all(f>=0 & f<=1 & theta>=0 & theta<1-1e-8) && all(Ts>=240 & Ts<=310), ...
    'ThreePhase:Domain','Phase fractions or solid temperatures exceed verification domain.');
mechanics=solve_cohesive_ring_mechanics(pg,Ts,model,p,history,loadingAtHistory);
u=mechanics.u;
gap=model.gapOperator*u; b=p.aperture0_m+p.wetAreaFraction*gap;
assert(all(b>0),'ThreePhase:Domain','Nonpositive available aperture.');
q=joint_threephase_gibbs(pg,T,p);
vc=(1-f).*q.volume(:,1)+f.*q.volume(:,2);
mc=(1-theta).*b./vc; mv=theta.*b./q.volume(:,3);
masses=[(1-f).*mc,f.*mc,mv]; volumes=masses.*q.volume;
energy=sum(masses.*q.energy,2); entropy=sum(masses.*q.entropy,2);
[~,stable]=min(q.g(:,1:2),[],2); ids=sub2ind(size(q.g),(1:n)',stable);
s=struct('x',[gap,sum(masses,2),masses(:,2),T],'u',u, ...
    'solidTemperature_K',Ts,'pressure',pg,'absolutePressure_Pa',pa, ...
    'fraction',f,'vaporFraction',theta,'phaseMass_kgm2',masses, ...
    'liquidVolume_m',volumes(:,1),'iceVolume_m',volumes(:,2),'vaporVolume_m',volumes(:,3), ...
    'internalEnergy_Jm2',energy,'physicalInternalEnergy_Jm2',energy-p.atmosphericPressure_Pa*b, ...
    'entropy_Jm2K',entropy,'freeEnergy_Jm2',sum(masses.*q.g,2)-pg.*b, ...
    'phaseAffinity_Jkg',q.g(:,2)-q.g(:,1),'vaporAffinity_Jkg',q.g(:,3)-q.g(ids), ...
    'aperture_m',b,'gibbs',q,'cohesive',mechanics.local);
if nargout<2, return; end
dp=zeros(n,nz); dp(:,il)=diag(pa);
df=zeros(n,nz); df(:,iff)=eye(n); dvf=zeros(n,nz); dvf(:,iv)=eye(n);
dT=zeros(n,nz); dT(:,it)=eye(n);
du=mechanics.pressureDerivative*dp; du(:,its)=mechanics.temperatureDerivative;
dg=model.gapOperator*du; db=p.wetAreaFraction*dg;
dvc=((1-f).*q.volumePressure(:,1)+f.*q.volumePressure(:,2)).*dp+ ...
    (q.volume(:,2)-q.volume(:,1)).*df;
dmc=((1-theta).*db-b.*dvf-mc.*dvc)./vc;
dmv=(theta.*db+b.*dvf-mv.*(q.volumePressure(:,3).*dp+q.volumeTemperature(:,3).*dT))./q.volume(:,3);
dml=(1-f).*dmc-mc.*df; dmi=f.*dmc+mc.*df;
dl=q.volume(:,1).*dml+masses(:,1).*q.volumePressure(:,1).*dp;
di=q.volume(:,2).*dmi+masses(:,2).*q.volumePressure(:,2).*dp;
de=q.energy(:,1).*dml+q.energy(:,2).*dmi+q.energy(:,3).*dmv+ ...
    sum(masses.*q.energyPressure,2).*dp+sum(masses.*q.energyTemperature,2).*dT;
dA=(q.volume(:,2)-q.volume(:,1)).*dp-(q.entropy(:,2)-q.entropy(:,1)).*dT;
dV=(q.volume(:,3)-q.volume(ids)).*dp-(q.entropy(:,3)-q.entropy(ids)).*dT;
d=struct('u',du,'gap',dg,'pressure',dp,'fraction',df,'vaporFraction',dvf, ...
    'temperature',dT,'mass',dmc+dmv,'ice',dmi,'vaporMass',dmv, ...
    'liquidVolume',dl,'iceVolume',di,'vaporVolume',db-dl-di, ...
    'energy',de,'affinity',dA,'vaporAffinity',dV);
end
