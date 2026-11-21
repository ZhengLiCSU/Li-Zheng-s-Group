function [s,d] = premelting_slit_state(pl,T,b,p,f)
%PREMELTING_SLIT_STATE Two wet wall films and one stationary ice slab.
% A=ml*g_l+mi*g_i-pl*Vl-pi*Vi+2*W(h), S=sum(m*s), U=A+T*S.
% At equilibrium dA=-S*dT-pi*db+mu*dM. Gauge energy includes ambient work.
% This is a pre-existing ice topology; it deliberately does not patch in
% liquid/vapour states, nucleation or disappearance of the two interfaces.
arguments
    pl (:,1) double {mustBeFinite}
    T (:,1) double {mustBeFinite}
    b (:,1) double {mustBeFinite,mustBePositive}
    p (1,1) struct = joint_threephase_parameters()
    f (1,1) struct = premelting_film_parameters()
end
assert(p.wetAreaFraction==1,'Premelting:Area','This resolved planar slit requires unit wetted area fraction.');
e=premelting_phase_equilibrium(pl,T,p,f); h=e.film_m;
assert(numel(b)==numel(pl),'Premelting:Shape','Aperture size mismatch.');
assert(all(2*h./b<f.maximumFilmFraction),'Premelting:IceTopology', ...
    'Opposing films interact or the ice core disappears; thin separated-film model invalid.');
vl=e.liquid.volume(:,1); vi=e.ice.volume(:,2);
Vl=2*h; Vi=b-Vl; ml=Vl./vl; mi=Vi./vi;
S=ml.*e.liquid.entropy(:,1)+mi.*e.ice.entropy(:,2);
U=ml.*e.liquid.energy(:,1)+mi.*e.ice.energy(:,2)+2*e.potential.energy_Jm2;
A=ml.*e.liquid.g(:,1)+mi.*e.ice.g(:,2)-pl.*Vl-e.icePressure_Pa.*Vi+2*e.potential.energy_Jm2;
% Two stationary no-slip films in parallel. Cubing TOTAL liquid thickness
% as one unobstructed slot would overestimate this conductance by four.
conductance=2*h.^3/(12*p.viscosity_Pas);
s=struct('equilibrium',e,'mass_kgm2',ml+mi,'phaseMass_kgm2',[ml,mi], ...
    'phaseVolume_m',[Vl,Vi],'aperture_m',b,'freeEnergy_Jm2',A, ...
    'internalEnergy_Jm2',U,'physicalInternalEnergy_Jm2',U-p.atmosphericPressure_Pa*b, ...
    'entropy_Jm2K',S,'filmEnergy_Jm2',2*e.potential.energy_Jm2, ...
    'volumeConductance_m3Pas',conductance,'chemicalMobility',conductance./vl.^2);
if nargout<2, return; end
% Derivative columns [liquid pressure Pa, temperature K, aperture m].
dh=[e.filmDerivative,zeros(size(pl))]; dpi=[e.icePressureDerivative,zeros(size(pl))];
dpl=repmat([1,0,0],numel(pl),1); dT=repmat([0,1,0],numel(pl),1); db=repmat([0,0,1],numel(pl),1);
dvl=e.liquid.volumePressure(:,1).*dpl;
dvi=e.ice.volumePressure(:,2).*dpi;
dml=(2*dh-ml.*dvl)./vl; dmi=(db-2*dh-mi.*dvi)./vi;
dU=e.liquid.energy(:,1).*dml+e.ice.energy(:,2).*dmi+ ...
    ml.*(e.liquid.energyPressure(:,1).*dpl+e.liquid.energyTemperature(:,1).*dT)+ ...
    mi.*(e.ice.energyPressure(:,2).*dpi+e.ice.energyTemperature(:,2).*dT)-2*e.pressureSplit_Pa.*dh;
dS=e.liquid.entropy(:,1).*dml+e.ice.entropy(:,2).*dmi+ ...
    (ml*p.heatWater_JkgK+mi*p.heatIce_JkgK)./T.*dT;
d=struct('mass',dml+dmi,'energy',dU,'entropy',dS, ...
    'freeEnergy',dU-T.*dS-S.*dT,'film',dh,'icePressure',dpi, ...
    'volumeConductance',h.^2/(2*p.viscosity_Pas).*dh);
end
