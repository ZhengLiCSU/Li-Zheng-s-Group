function [x,state] = joint_closed_equilibrium(mass,temperature,p)
%JOINT_CLOSED_EQUILIBRIUM Exact local phase minimum at fixed total mass.
% A scalar bracketed phase solve plus normal equilibrium. Used as a Newton
% predictor only: an open interface still has to solve its complete mass PDE.
arguments
    mass (:,1) double {mustBePositive,mustBeFinite}
    temperature (1,1) double {mustBePositive,mustBeFinite}
    p (1,1) struct = joint_thermodynamic_parameters()
end
n = numel(mass); x = zeros(n,4); primaryPressure = zeros(n,1);
for k = 1:n
    lower = affinity(0,mass(k),temperature,p);
    upper = affinity(1,mass(k),temperature,p);
    if lower>=0
        fraction = 0;
    elseif upper<=0
        fraction = 1;
    else
        fraction = fzero(@(f)affinity(f,mass(k),temperature,p),[0,1]);
    end
    [pressure,gap] = mechanicalPressure(fraction,mass(k),p);
    x(k,:) = [gap,mass(k),fraction*mass(k),temperature];
    primaryPressure(k) = pressure;
end
state = joint_helmholtz_state(x,p,true,primaryPressure);
end

function value = affinity(fraction,mass,temperature,p)
pressure = mechanicalPressure(fraction,mass,p);
shape = temperature-p.meltingTemperature_K-temperature*log(temperature/p.meltingTemperature_K);
value = (p.heatIce_JkgK-p.heatWater_JkgK)*shape+ ...
    p.latentHeat_Jkg*(temperature/p.meltingTemperature_K-1)+p.capillaryShift_Jkg+ ...
    (1/p.rhoIce-1/p.rhoWater)*pressure-0.5*( ...
    1/(p.rhoIce*p.bulkIce_Pa)-1/(p.rhoWater*p.bulkWater_Pa))*pressure^2;
end

function [pressure,gap] = mechanicalPressure(fraction,mass,p)
v = (1-fraction)/p.rhoWater+fraction/p.rhoIce;
c = (1-fraction)/(p.rhoWater*p.bulkWater_Pa)+fraction/(p.rhoIce*p.bulkIce_Pa);
stiffness = p.normalStiffness_Pam;
preload = p.normalStiffness_Pam*p.referenceGap_m+p.externalTraction_Pa;
pressure = (mass*v-p.aperture0_m-p.wetAreaFraction*preload/stiffness)/ ...
    (mass*c+p.wetAreaFraction^2/stiffness);
gap = (p.wetAreaFraction*pressure+preload)/stiffness;
if gap<0
    stiffness = stiffness+p.contactStiffness_Pam;
    pressure = (mass*v-p.aperture0_m-p.wetAreaFraction*preload/stiffness)/ ...
        (mass*c+p.wetAreaFraction^2/stiffness);
    gap = (p.wetAreaFraction*pressure+preload)/stiffness;
end
end
