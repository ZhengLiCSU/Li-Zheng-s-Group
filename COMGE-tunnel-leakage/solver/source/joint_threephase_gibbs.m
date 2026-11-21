function q = joint_threephase_gibbs(pressure,temperature,p)
%JOINT_THREEPHASE_GIBBS One gauge-referenced potential for all three phases.
% Columns: liquid, ice, vapour. g_p=v, -g_T=s, u=g+Ts-p_g v.
% Physical internal energy equals u-p_atm*v. Using u includes the constant
% ambient-pressure work reservoir, consistent with gauge-pressure mechanics.
% No pressure floor or artificial refill is used.
arguments
    pressure (:,1) double {mustBeFinite}
    temperature (:,1) double {mustBeFinite}
    p (1,1) struct = joint_threephase_parameters()
end
assert(numel(pressure)==numel(temperature),'ThreePhase:Shape','State sizes must match.');
pa=pressure+p.atmosphericPressure_Pa; T=temperature; pg=pressure;
assert(all(pa>0 & T>=240 & T<=310),'ThreePhase:Domain', ...
    'Require positive absolute pressure and 240<=T<=310 K for this local EOS.');
assert(all(abs(pg)<=p.maximumCompression*min(p.bulkWater_Pa,p.bulkIce_Pa)), ...
    'ThreePhase:Domain','Condensed small-compression domain exceeded.');
v0=[1/p.rhoWater,1/p.rhoIce];
c=v0./[p.bulkWater_Pa,p.bulkIce_Pa];
cp=[p.heatWater_JkgK,p.heatIce_JkgK,p.vaporHeatCapacity_JkgK];
Tm=p.meltingTemperature_K; R=p.vaporGasConstant_JkgK; L=p.vaporLatentHeat_Jkg;
logT=log(T/Tm); shape=T-Tm-T.*logT;
g=zeros(numel(T),3); entropy=g; v=g; vp=g; vT=g;
g(:,1:2)=shape.*cp(1:2)+pg.*v0-0.5*pg.^2.*c;
g(:,2)=g(:,2)+p.latentHeat_Jkg*(T/Tm-1)+p.iceReferenceShift_Jkg;
g(:,3)=p.vaporReference_Jkg+cp(3)*shape+L*(1-T/Tm)+R*T.*log(pa/p.triplePressure_Pa);
entropy(:,1:2)=logT.*cp(1:2);
entropy(:,2)=entropy(:,2)-p.latentHeat_Jkg/Tm;
entropy(:,3)=cp(3)*logT+L/Tm-R*log(pa/p.triplePressure_Pa);
v(:,1:2)=v0-pg.*c; v(:,3)=R*T./pa;
vp(:,1:2)=-repmat(c,numel(T),1); vp(:,3)=-R*T./pa.^2;
vT(:,3)=R./pa;
u=g+T.*entropy-pg.*v;
up=-T.*vT-pg.*vp; uT=cp-pg.*vT;
q=struct('g',g,'entropy',entropy,'volume',v,'volumePressure',vp, ...
    'volumeTemperature',vT,'energy',u,'energyPressure',up,'energyTemperature',uT);
end
