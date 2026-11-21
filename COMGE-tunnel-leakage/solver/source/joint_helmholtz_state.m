function s = joint_helmholtz_state(x,p,enforceDomain,primaryPressure)
%JOINT_HELMHOLTZ_STATE Common Helmholtz energy and exact local tangent.
% x(:,1:4) = [gap increment m, total mass kg/m2, ice mass kg/m2, T K].
% All extensive quantities are per reference wall area (unit axial width).
% Liquid transport changes total mass at fixed ice mass. Its conjugate is
% mu_l = dPsi/dM|I. Pressure and normal traction share the SAME aperture.
% Quadratic compressibility is a small-compression approximation, with an
% explicit domain guard. It does not represent cavitation or ice lenses.
arguments
    x (:,4) double {mustBeReal,mustBeFinite}
    p (1,1) struct = joint_thermodynamic_parameters()
    enforceDomain (1,1) logical = true
    primaryPressure (:,1) double {mustBeReal,mustBeFinite} = zeros(0,1)
end
gap = x(:,1); mass = x(:,2); ice = x(:,3); temperature = x(:,4);
assert(all(mass>0 & ice>=0 & ice<=mass & temperature>0), ...
    'JointThermo:StateDomain','Require M>0, 0<=I<=M and absolute T>0.');
vw = 1/p.rhoWater; vi = 1/p.rhoIce;
cw = vw/p.bulkWater_Pa; ci = vi/p.bulkIce_Pa;
liquid = mass-ice;
volume = p.aperture0_m+p.wetAreaFraction*gap;
assert(all(volume>0),'JointThermo:ApertureDomain','Available aperture must be positive.');
naturalVolume = liquid*vw+ice*vi;
compliance = liquid*cw+ice*ci;
excessVolume = (mass-p.rhoWater*p.aperture0_m)*vw+ice*(vi-vw)-p.wetAreaFraction*gap;
pressure = excessVolume./compliance;
if ~isempty(primaryPressure)
    assert(numel(primaryPressure)==size(x,1),'JointThermo:PressureShape','One pressure per cell is required.');
    volumeError = naturalVolume-volume-compliance.*primaryPressure;
    assert(all(abs(volumeError)<=64*eps*(naturalVolume+volume)), ...
        'JointThermo:PressureClosure','Primary pressure must satisfy the same volume constraint.');
    pressure = primaryPressure;
end
liquidVolume = liquid.*(vw-cw*pressure);
iceVolume = ice.*(vi-ci*pressure);
if enforceDomain
    assert(all(pressure>=p.minimumPressure_Pa),'JointThermo:CavitationDomain', ...
        'Saturated single-pressure branch cannot represent this tensile pressure.');
    assert(all(abs(pressure)/min(p.bulkWater_Pa,p.bulkIce_Pa)<=p.maximumCompression), ...
        'JointThermo:CompressionDomain','Small-compression domain exceeded.');
end
logTemperature = log(temperature/p.meltingTemperature_K);
thermalShape = temperature-p.meltingTemperature_K-temperature.*logTemperature;
fLiquid = p.heatWater_JkgK*thermalShape;
fIce = p.heatIce_JkgK*thermalShape+ ...
    p.latentHeat_Jkg*(temperature/p.meltingTemperature_K-1)+p.capillaryShift_Jkg;
heatCapacity = liquid*p.heatWater_JkgK+ice*p.heatIce_JkgK;
compressionEnergy = 0.5*compliance.*pressure.^2;
contactGap = min(gap,0);
mechanicalEnergy = 0.5*p.normalStiffness_Pam*(gap-p.referenceGap_m).^2+ ...
    0.5*p.contactStiffness_Pam*contactGap.^2;
externalPotential = -p.externalTraction_Pa*gap;
freeEnergy = liquid.*fLiquid+ice.*fIce+compressionEnergy+ ...
    mechanicalEnergy+externalPotential;
entropy = heatCapacity.*logTemperature-ice*p.latentHeat_Jkg/p.meltingTemperature_K;
internalEnergy = heatCapacity.*(temperature-p.meltingTemperature_K)- ...
    ice*p.latentHeat_Jkg+ice*p.capillaryShift_Jkg+compressionEnergy+mechanicalEnergy;
volumeGradient = [-p.wetAreaFraction,vw,vi-vw];
complianceGradient = [0,cw,ci-cw];
a = volumeGradient-pressure.*complianceGradient;
pressureGradient = a./compliance;
gradient = zeros(size(x));
gradient(:,1:3) = pressure.*volumeGradient-0.5*pressure.^2.*complianceGradient;
gradient(:,1) = gradient(:,1)+p.normalStiffness_Pam*(gap-p.referenceGap_m)+ ...
    p.contactStiffness_Pam*contactGap-p.externalTraction_Pa;
gradient(:,2) = gradient(:,2)+fLiquid;
gradient(:,3) = gradient(:,3)+fIce-fLiquid;
gradient(:,4) = -entropy;
n = size(x,1);
hessian = zeros(4,4,n);
for k = 1:n
    hessian(1:3,1:3,k) = (a(k,:)'*a(k,:))/compliance(k);
    hessian(1,1,k) = hessian(1,1,k)+p.normalStiffness_Pam+ ...
        p.contactStiffness_Pam*(gap(k)<0);
    hessian(2,4,k) = -p.heatWater_JkgK*logTemperature(k);
    hessian(3,4,k) = -(p.heatIce_JkgK-p.heatWater_JkgK)*logTemperature(k)+ ...
        p.latentHeat_Jkg/p.meltingTemperature_K;
    hessian(4,2:3,k) = hessian(2:3,4,k)';
    hessian(4,4,k) = -heatCapacity(k)/temperature(k);
end
liquidVolumeGradient = -liquid.*cw.*pressureGradient;
liquidVolumeGradient(:,2) = liquidVolumeGradient(:,2)+vw-cw*pressure;
liquidVolumeGradient(:,3) = liquidVolumeGradient(:,3)-vw+cw*pressure;
energyGradient = gradient;
energyGradient(:,1) = energyGradient(:,1)+p.externalTraction_Pa;
energyGradient(:,2) = energyGradient(:,2)+temperature.*p.heatWater_JkgK.*logTemperature;
energyGradient(:,3) = energyGradient(:,3)+temperature.*( ...
    (p.heatIce_JkgK-p.heatWater_JkgK)*logTemperature-p.latentHeat_Jkg/p.meltingTemperature_K);
energyGradient(:,4) = heatCapacity;
s = struct('freeEnergy_Jm2',freeEnergy,'internalEnergy_Jm2',internalEnergy, ...
    'entropy_Jm2K',entropy,'heatCapacity_Jm2K',heatCapacity, ...
    'gradient',gradient,'hessian',hessian,'pressure_Pa',pressure, ...
    'pressureGradient',pressureGradient,'liquidChemicalPotential_Jkg',gradient(:,2), ...
    'phaseAffinity_Jkg',gradient(:,3),'normalResidual_Pa',gradient(:,1), ...
    'fluidNormalTraction_Pa',-p.wetAreaFraction*pressure, ...
    'aperture_m',volume,'liquidVolume_m',liquidVolume,'iceVolume_m',iceVolume, ...
    'liquidVolumeGradient',liquidVolumeGradient,'energyGradient',energyGradient, ...
    'compressionEnergy_Jm2',compressionEnergy,'mechanicalEnergy_Jm2',mechanicalEnergy, ...
    'liquidMass_kgm2',liquid,'iceMass_kgm2',ice);
end
