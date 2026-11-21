function [residual,tangent,flow,state] = joint_interface_residual(x,old,temperature,dt,settings,p,primaryPressure)
%JOINT_INTERFACE_RESIDUAL Monolithic mass, phase KKT and normal equilibrium.
% Cell-centred conservative transport along a planar slit. All cells share
% a prescribed thermostat temperature; this is NOT a full heat PDE solver.
% The Onsager mobility rho_l^2*b_l^3/(12*mu) gives planar Poiseuille flow in
% the continuum limit. No fitted damage multiplier is used. Ice is immobile.
arguments
    x (:,3) double {mustBeReal,mustBeFinite}
    old (:,4) double {mustBeReal,mustBeFinite}
    temperature (1,1) double {mustBePositive,mustBeFinite}
    dt (1,1) double {mustBePositive,mustBeFinite}
    settings (1,1) struct
    p (1,1) struct
    primaryPressure (:,1) double {mustBeReal,mustBeFinite} = zeros(0,1)
end
n = size(x,1);
assert(isequal(size(old),[n,4]),'JointInterface:Shape','Old state size must match cells.');
state = joint_helmholtz_state([x,repmat(temperature,n,1)],p,false,primaryPressure);
dx = settings.length_m/n;
massScale = p.rhoWater*p.aperture0_m;
gapScale = p.aperture0_m/p.wetAreaFraction;
stressScale = p.normalStiffness_Pam*gapScale;
affinityScale = p.latentHeat_Jkg/p.meltingTemperature_K;
indices = reshape(1:3*n,3,n)';
mu = state.liquidChemicalPotential_Jkg;
muGradient = zeros(n,3);
for k = 1:n
    muGradient(k,:) = state.hessian(2,1:3,k);
end
vw = 1/p.rhoWater; cw = vw/p.bulkWater_Pa;
liquidSpecificVolume = vw-cw*state.pressure_Pa;
% Volume is per TOTAL reference area. The physical liquid-layer thickness
% on the wet part is b_l/chi, so chi*(b_l/chi)^3 = b_l^3/chi^2.
mobility = state.liquidVolume_m.^3./(12*p.viscosity_Pas*liquidSpecificVolume.^2*p.wetAreaFraction^2);
mobilityGradient = state.liquidVolume_m.^2./(4*p.viscosity_Pas*liquidSpecificVolume.^2*p.wetAreaFraction^2).* ...
    state.liquidVolumeGradient+state.liquidVolume_m.^3*cw./ ...
    (6*p.viscosity_Pas*liquidSpecificVolume.^3*p.wetAreaFraction^2).*state.pressureGradient;
flux = zeros(n+1,1); fluxTangent = sparse(n+1,3*n);
potentialDifference = zeros(n+1,1);
reservoirMu = reservoirPotential(settings.reservoirPressure_Pa,temperature,p);
if ~settings.sealedEnds(1)
    potentialDifference(1) = pressurePotentialDifference( ...
        settings.reservoirPressure_Pa(1),state.pressure_Pa(1),vw,cw);
    conductance = 2*mobility(1)/dx;
    flux(1) = conductance*potentialDifference(1);
    fluxTangent(1,indices(1,:)) = 2*mobilityGradient(1,:)/dx*potentialDifference(1)- ...
        conductance*muGradient(1,:);
end
for face = 2:n
    left = face-1; right = face;
    denominator = mobility(left)+mobility(right);
    potentialDifference(face) = pressurePotentialDifference( ...
        state.pressure_Pa(left),state.pressure_Pa(right),vw,cw);
    if denominator>0
        conductance = 2*mobility(left)*mobility(right)/(dx*denominator);
        derivativeLeft = 2*mobility(right)^2/(dx*denominator^2)*mobilityGradient(left,:);
        derivativeRight = 2*mobility(left)^2/(dx*denominator^2)*mobilityGradient(right,:);
        flux(face) = conductance*potentialDifference(face);
        fluxTangent(face,indices(left,:)) = derivativeLeft*potentialDifference(face)+ ...
            conductance*muGradient(left,:);
        fluxTangent(face,indices(right,:)) = derivativeRight*potentialDifference(face)- ...
            conductance*muGradient(right,:);
    end
end
if ~settings.sealedEnds(2)
    potentialDifference(end) = pressurePotentialDifference( ...
        state.pressure_Pa(end),settings.reservoirPressure_Pa(2),vw,cw);
    conductance = 2*mobility(end)/dx;
    flux(end) = conductance*potentialDifference(end);
    fluxTangent(end,indices(end,:)) = 2*mobilityGradient(end,:)/dx*potentialDifference(end)+ ...
        conductance*muGradient(end,:);
    if settings.outletSeepage && potentialDifference(end)<=0
        flux(end) = 0;
        fluxTangent(end,:) = 0;
    end
end
residual = zeros(3*n,1); tangent = sparse(3*n,3*n);
for k = 1:n
    ids = indices(k,:);
    residual(ids(1)) = state.normalResidual_Pa(k)/stressScale;
    tangent(ids(1),ids) = state.hessian(1,1:3,k)/stressScale;
    residual(ids(2)) = (x(k,2)-old(k,2)+dt/dx*(flux(k+1)-flux(k)))/massScale;
    tangent(ids(2),:) = dt/(dx*massScale)*(fluxTangent(k+1,:)-fluxTangent(k,:));
    tangent(ids(2),ids(2)) = tangent(ids(2),ids(2))+1/massScale;
    affinity = state.phaseAffinity_Jkg(k)+p.phaseViscosity_Jm2skg2*(x(k,3)-old(k,3))/dt;
    affinityGradient = state.hessian(3,1:3,k);
    affinityGradient(3) = affinityGradient(3)+p.phaseViscosity_Jm2skg2/dt;
    trialIce = x(k,3)/massScale-affinity/affinityScale;
    if trialIce<=0
        residual(ids(3)) = x(k,3)/massScale;
        tangent(ids(3),ids(3)) = 1/massScale;
    elseif trialIce>=x(k,2)/massScale
        residual(ids(3)) = (x(k,3)-x(k,2))/massScale;
        tangent(ids(3),ids(2:3)) = [-1,1]/massScale;
    else
        residual(ids(3)) = affinity/affinityScale;
        tangent(ids(3),ids) = affinityGradient/affinityScale;
    end
end
flow = struct('massFlux_kgms',flux,'referenceVolumeFlux_m2s',flux/p.rhoWater, ...
    'potentialDifference_Jkg',potentialDifference,'reservoirPotential_Jkg',reservoirMu, ...
    'cellMobility',mobility,'fluxTangent',fluxTangent);
end

function difference = pressurePotentialDifference(left,right,vw,cw)
% Uniform T: the thermal reference cancels analytically, not numerically.
% Avoid subtracting nearly equal chemical potentials with a large common
% thermal offset. This preserves the same thermodynamic driving force.
difference = (vw-0.5*cw*(left+right))*(left-right);
end

function chemicalPotential = reservoirPotential(pressure,temperature,p)
vw = 1/p.rhoWater; cw = vw/p.bulkWater_Pa;
fLiquid = p.heatWater_JkgK*(temperature-p.meltingTemperature_K- ...
    temperature*log(temperature/p.meltingTemperature_K));
chemicalPotential = fLiquid+vw*pressure-0.5*cw*pressure.^2;
end
