function [next,report] = joint_adiabatic_equilibrium(initial,p)
%JOINT_ADIABATIC_EQUILIBRIUM Sealed cell with unknown equilibrium temperature.
% Conserves the total internal energy (fluid + elastic wall) and total mass.
% This is an equilibrium recalescence benchmark, not a finite-rate heat PDE.
arguments
    initial (1,4) double {mustBeReal,mustBeFinite}
    p (1,1) struct = joint_thermodynamic_parameters()
end
assert(p.externalTraction_Pa==0,'JointThermo:AdiabaticLoading', ...
    'This closed-energy benchmark requires zero external traction.');
before = joint_helmholtz_state(initial,p);
mass = initial(2); targetEnergy = before.internalEnergy_Jm2;
bracket = [min(initial(4),p.meltingTemperature_K)-20, ...
    max(initial(4),p.meltingTemperature_K)+20];
temperature = fzero(@(t)energyResidual(t,mass,targetEnergy,p),bracket, ...
    optimset('TolX',eps));
[next,after] = joint_closed_equilibrium(mass,temperature,p);
energyDefect = after.internalEnergy_Jm2-targetEnergy;
entropyProduction = after.entropy_Jm2K-before.entropy_Jm2K;
assert(abs(energyDefect)<=1e-8*max(1,abs(targetEnergy)) && entropyProduction>=-1e-9, ...
    'JointThermo:AdiabaticAudit','Energy defect %.4g J/m2, entropy production %.4g J/(m2 K).', ...
    energyDefect,entropyProduction);
report = struct('energyDefect_Jm2',energyDefect,'massDefect_kgm2',next(2)-mass, ...
    'entropyProduction_Jm2K',entropyProduction,'initial',before,'final',after);
end

function value = energyResidual(temperature,mass,targetEnergy,p)
[~,state] = joint_closed_equilibrium(mass,temperature,p);
value = state.internalEnergy_Jm2-targetEnergy;
end
