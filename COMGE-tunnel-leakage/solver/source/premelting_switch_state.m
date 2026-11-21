function [s,d] = premelting_switch_state(pl,T,b,theta,p,f)
%PREMELTING_SWITCH_STATE Conserved water in liquid / wet-ice / vapour areas.
% theta contains AREA fractions, not phase mass or volume fractions.
% Omega_i=-pi*b+2*(pi-pl)*h+2*(gamma_iw+W(h)). The film minimizes Omega_i.
% F=sum(m*g-p*V)+surface, U=F+T*S; no reset when a phase disappears.
% Derivative columns [pl,T,b,theta_l,theta_i,theta_v]. Fractions may leave
% the simplex during algebraic Newton trials; accepted states must pass KKT.
arguments
    pl (:,1) double {mustBeFinite}
    T (:,1) double {mustBeFinite}
    b (:,1) double {mustBeFinite,mustBePositive}
    theta (:,3) double {mustBeFinite}
    p (1,1) struct = joint_threephase_parameters()
    f (1,1) struct = premelting_switch_parameters()
end
n=numel(pl); assert(numel(b)==n && size(theta,1)==n && p.wetAreaFraction==1,'PhaseSwitch:Shape','Invalid shape or wet area.');
assert(all(b>2*f.minimumFilm_m),'PhaseSwitch:Domain','Aperture below molecular film domain.');
q=premelting_switch_eos(pl,T,p); split=q.pressure(:,2)-pl;
h=zeros(n,1); branch=zeros(n,1);
lo=f.minimumFilm_m; hi=b/2; stableHi=min(hi,f.spinodalFilm_m);
vl=premelting_film_potential(lo,f); vh=premelting_film_potential(stableHi,f);
rootMask=split>vh.pressure_Pa & split<vl.pressure_Pa;
roots=zeros(n,1);
if any(rootMask)
    roots(rootMask)=solve_premelting_film_roots(split(rootMask),lo,stableHi(rootMask),f);
end
% Evaluate the exact same lower, upper and stable interior candidates in
% one batch. Missing interior candidates have infinite comparison energy.
candidates=[lo*ones(n,1),hi,lo*ones(n,1)]; candidates(rootMask,3)=roots(rootMask);
potential=premelting_film_potential(candidates(:),f);
energies=2*split.*candidates+2*reshape(potential.energy_Jm2,n,3);
energies(~rootMask,3)=Inf;
[~,which]=min(energies,[],2);
h=candidates(sub2ind(size(candidates),(1:n)',which));
types=[-1;1;0]; branch=types(which);
v=premelting_film_potential(h,f); dh=zeros(n,6);
interior=branch==0; upper=branch==1;
dh(interior,1)=(q.pressureDerivative(interior,2,1)-1)./v.pressureDerivative_Pam(interior);
dh(interior,2)=q.pressureDerivative(interior,2,2)./v.pressureDerivative_Pam(interior);
dh(upper,3)=0.5;
columnP=q.pressure; columnP(upper,2)=pl(upper)+v.pressure_Pa(upper);
omega=[-pl.*b,-q.pressure(:,2).*b+2*split.*h+2*(f.iceWaterSurfaceEnergy_Jm2+v.energy_Jm2), ...
    -q.pressure(:,3).*b+2*f.wallVaporExcess_Jm2];
V=[theta(:,1).*b+2*theta(:,2).*h,theta(:,2).*(b-2*h),theta(:,3).*b];
mass=V./q.volume;
surface=2*theta(:,2).*(f.iceWaterSurfaceEnergy_Jm2+v.energy_Jm2)+2*theta(:,3)*f.wallVaporExcess_Jm2;
S=sum(mass.*q.entropy,2); U=sum(mass.*q.energy,2)+surface;
A=sum(mass.*q.g-q.pressure.*V,2)+surface;
viscosity=p.viscosity_Pas*ones(n,1); dviscosity=zeros(n,6);
if isfield(p,'viscosityModel')
    assert(ismember(p.viscosityModel,{'constant','iapws2008'}),'PhaseSwitch:Material','Unsupported water viscosity model.');
    if strcmp(p.viscosityModel,'iapws2008')
        [viscosity,vt,vr]=water_viscosity_iapws2008(T,1./q.volume(:,1));
        dviscosity(:,1)=-vr.*q.volumePressure(:,1)./q.volume(:,1).^2;
        dviscosity(:,2)=vt-vr.*q.volumeTemperature(:,1)./q.volume(:,1).^2;
    end
end
C=(theta(:,1).*b.^3+2*theta(:,2).*h.^3)./(12*viscosity);
s=struct('areaFraction',theta,'phaseMass_kgm2',mass,'phaseVolume_m',V,'mass_kgm2',sum(mass,2), ...
    'freeEnergy_Jm2',A,'internalEnergy_Jm2',U,'entropy_Jm2K',S,'surfaceEnergy_Jm2',surface, ...
    'physicalInternalEnergy_Jm2',U-p.atmosphericPressure_Pa*b,'aperture_m',b, ...
    'film_m',h,'filmBranch',branch,'wallPressure_Pa',sum(theta.*columnP,2),'grandPotential_Jm2',omega, ...
    'volumeConductance_m3Pas',C,'chemicalMobility',C./q.volume(:,1).^2,'eos',q, ...
    'waterViscosity_Pas',viscosity);
if nargout<2, return; end
db=repmat([0,0,1,0,0,0],n,1); dT=repmat([0,1,0,0,0,0],n,1);
dth=zeros(n,3,6); dth(:,1,4)=1; dth(:,2,5)=1; dth(:,3,6)=1;
dV=zeros(n,3,6); dm=dV; dp=dV; dOmega=dV; dColumnP=dV;
dp(:,:,1:2)=q.pressureDerivative;
dV(:,1,:)=reshape(theta(:,1).*db+b.*reshape(dth(:,1,:),n,6)+2*theta(:,2).*dh+2*h.*reshape(dth(:,2,:),n,6),n,1,6);
dV(:,2,:)=reshape(theta(:,2).*(db-2*dh)+(b-2*h).*reshape(dth(:,2,:),n,6),n,1,6);
dV(:,3,:)=reshape(theta(:,3).*db+b.*reshape(dth(:,3,:),n,6),n,1,6);
dS=zeros(n,6); dU=dS; dP=dS;
dSurface=2*(f.iceWaterSurfaceEnergy_Jm2+v.energy_Jm2).*reshape(dth(:,2,:),n,6)- ...
    2*theta(:,2).*v.pressure_Pa.*dh+2*f.wallVaporExcess_Jm2*reshape(dth(:,3,:),n,6);
for k=1:3
    dpk=reshape(dp(:,k,:),n,6);
    dvk=q.volumePressure(:,k).*dpk+q.volumeTemperature(:,k).*dT;
    dmk=(reshape(dV(:,k,:),n,6)-mass(:,k).*dvk)./q.volume(:,k);
    dm(:,k,:)=reshape(dmk,n,1,6);
    cp=[p.heatWater_JkgK,p.heatIce_JkgK,p.vaporHeatCapacity_JkgK];
    dS=dS+q.entropy(:,k).*dmk+mass(:,k).*(-q.volumeTemperature(:,k).*dpk+cp(k)./T.*dT);
    dU=dU+q.energy(:,k).*dmk+mass(:,k).*(q.energyPressure(:,k).*dpk+q.energyTemperature(:,k).*dT);
    dc=dpk;
    if k==2 && any(upper)
        dc(upper,:)=reshape(dp(upper,1,:),sum(upper),6)+v.pressureDerivative_Pam(upper).*dh(upper,:);
    end
    dColumnP(:,k,:)=reshape(dc,n,1,6);
    dP=dP+theta(:,k).*dc+columnP(:,k).*reshape(dth(:,k,:),n,6);
end
dOmega(:,1,:)=reshape(-b.*reshape(dp(:,1,:),n,6)-pl.*db,n,1,6);
dOmega(:,2,:)=reshape((2*h-b).*reshape(dp(:,2,:),n,6)-2*h.*reshape(dp(:,1,:),n,6)-columnP(:,2).*db,n,1,6);
dOmega(:,3,:)=reshape(-b.*reshape(dp(:,3,:),n,6)-q.pressure(:,3).*db,n,1,6);
dC=(b.^3.*reshape(dth(:,1,:),n,6)+3*theta(:,1).*b.^2.*db+ ...
    2*h.^3.*reshape(dth(:,2,:),n,6)+6*theta(:,2).*h.^2.*dh)./(12*viscosity)-C./viscosity.*dviscosity;
dU=dU+dSurface; dMass=reshape(sum(dm,2),n,6);
d=struct('mass',dMass,'energy',dU,'entropy',dS,'freeEnergy',dU-T.*dS-S.*dT, ...
    'wallPressure',dP,'grandPotential',dOmega,'film',dh,'phaseMass',dm,'phaseVolume',dV, ...
    'volumeConductance',dC,'chemicalMobility',dC./q.volume(:,1).^2- ...
    2*C.*q.volumePressure(:,1)./q.volume(:,1).^3.*reshape(dp(:,1,:),n,6));
end


