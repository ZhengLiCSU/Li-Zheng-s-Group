function q = premelting_switch_eos(pl,T,p)
%PREMELTING_SWITCH_EOS Common chemical potential, distinct phase pressures.
% The liquid may sustain tension. Only VAPOUR requires positive absolute
% pressure. All reference constants match joint_threephase_gibbs exactly.
arguments
    pl (:,1) double {mustBeFinite}
    T (:,1) double {mustBeFinite}
    p (1,1) struct = joint_threephase_parameters()
end
assert(numel(pl)==numel(T) && all(T>=240 & T<=310),'PhaseSwitch:Domain','Invalid temperature or size.');
v0=[1/p.rhoWater,1/p.rhoIce]; c=v0./[p.bulkWater_Pa,p.bulkIce_Pa];
cp=[p.heatWater_JkgK,p.heatIce_JkgK,p.vaporHeatCapacity_JkgK];
Tm=p.meltingTemperature_K; R=p.vaporGasConstant_JkgK; L=p.vaporLatentHeat_Jkg;
lt=log(T/Tm); shape=T-Tm-T.*lt;
mu=cp(1)*shape+v0(1)*pl-0.5*c(1)*pl.^2;
iceBase=cp(2)*shape+p.latentHeat_Jkg*(T/Tm-1)+p.iceReferenceShift_Jkg;
D=mu-iceBase; discriminant=v0(2)^2-2*c(2)*D;
assert(all(discriminant>0),'PhaseSwitch:Domain','Ice EOS inversion exceeded its domain.');
pi=2*D./(v0(2)+sqrt(discriminant));
logPv=(mu-p.vaporReference_Jkg-cp(3)*shape-L*(1-T/Tm))./(R*T);
pa=p.triplePressure_Pa*exp(logPv); pv=pa-p.atmosphericPressure_Pa;
pressure=[pl,pi,pv];
assert(all(abs(pressure(:,1:2))<=p.maximumCompression*min(p.bulkWater_Pa,p.bulkIce_Pa),'all') && all(pa>0 & isfinite(pa)), ...
    'PhaseSwitch:Domain','Condensed EOS strain or vapour pressure domain exceeded.');
v=[v0(1)-c(1)*pl,v0(2)-c(2)*pi,R*T./pa];
vp=[-c(1)*ones(size(T)),-c(2)*ones(size(T)),-R*T./pa.^2];
vT=[zeros(numel(T),2),R./pa];
s=[cp(1)*lt,cp(2)*lt-p.latentHeat_Jkg/Tm,cp(3)*lt+L/Tm-R*logPv];
g=[mu,iceBase+v0(2)*pi-0.5*c(2)*pi.^2, ...
    p.vaporReference_Jkg+cp(3)*shape+L*(1-T/Tm)+R*T.*logPv];
u=g+T.*s-pressure.*v;
dp=zeros(numel(T),3,2); dp(:,:,1)=v(:,1)./v;
dp(:,:,2)=(s-s(:,1))./v;
q=struct('pressure',pressure,'pressureDerivative',dp,'g',g,'chemicalPotential_Jkg',mu, ...
    'entropy',s,'volume',v,'volumePressure',vp,'volumeTemperature',vT, ...
    'energy',u,'energyPressure',-T.*vT-pressure.*vp,'energyTemperature',cp-pressure.*vT);
end
