function p = joint_threephase_parameters()
%JOINT_THREEPHASE_PARAMETERS Pure-water equilibrium verification material.
% Separate from the saturated v16 model. Planar phases, no air, nucleation
% barrier, surface energy, ice shear stress, or vapour transport are claimed.
% Constant heat capacities / linear condensed EOS are local approximations.
p=joint_thermodynamic_parameters(struct('meltingTemperature_K',273.16));
p.atmosphericPressure_Pa=101325;
p.triplePressure_Pa=611.657;
p.vaporGasConstant_JkgK=461.51805;
p.vaporHeatCapacity_JkgK=1859;
p.vaporLatentHeat_Jkg=2.5009e6;
p.vaporConductivity_WmK=0.018;
% Calibrate chemical potentials, not a fitted cavitation threshold.
% The liquid reference is unchanged. Ice and vapour meet it at the IAPWS
% triple point. This offset is a reference correction, NOT pore capillarity.
pr=p.triplePressure_Pa-p.atmosphericPressure_Pa;
vw=1/p.rhoWater; vi=1/p.rhoIce;
cw=vw/p.bulkWater_Pa; ci=vi/p.bulkIce_Pa;
p.iceReferenceShift_Jkg=-(vi-vw)*pr+0.5*(ci-cw)*pr^2;
p.vaporReference_Jkg=vw*pr-0.5*cw*pr^2;
end
