function [flux,derivative,audit] = joint_nonisothermal_face(left,right,mobility,thermal,p)
%JOINT_NONISOTHERMAL_FACE Conservative, entropy-producing liquid/heat port.
% left/right = [pressure Pa, T K]. Mobility already includes face distance;
% thermal is Fourier conductance W/(m K). Ice is stationary. Derivative
% columns are [pL,TL,pR,TR,mobility,thermal]. This neglects thermodiffusion
% and gravity; a temperature gradient at equal pressure drives heat only.
arguments
    left (1,2) double {mustBeFinite}
    right (1,2) double {mustBeFinite}
    mobility (1,1) double {mustBeNonnegative,mustBeFinite}
    thermal (1,1) double {mustBeNonnegative,mustBeFinite}
    p (1,1) struct = joint_thermodynamic_parameters()
end
assert(left(2)>0 && right(2)>0,'JointNonisothermal:Temperature','Absolute temperatures must be positive.');
vw = 1/p.rhoWater; cw = vw/p.bulkWater_Pa;
pressure = [left(1),right(1)]; temperature = [left(2),right(2)];
P = vw*pressure-0.5*cw*pressure.^2;
dP = vw-cw*pressure;
deltaP = (vw-0.5*cw*sum(pressure))*(pressure(1)-pressure(2));
deltaT = temperature(1)-temperature(2);
mass = mobility*deltaP;
% Upwind sensible enthalpy plus centred pressure enthalpy. This admits an
% EXACT nonnegative face-entropy decomposition, including at unequal T.
up = 1+(mass<0);
h = p.heatWater_JkgK*(temperature(up)-p.meltingTemperature_K)+mean(P);
energy = mass*h+thermal*deltaT;
flux = [mass;energy];
if nargout<2, return; end
dmass = [mobility*dP(1),0,-mobility*dP(2),0,deltaP,0];
dh = [0.5*dP(1),0,0.5*dP(2),0,0,0]; dh(2*up) = p.heatWater_JkgK;
derivative = [dmass;h*dmass+mass*dh+[0,thermal,0,-thermal,0,deltaT]];
if nargout<3, return; end
% Inverse logarithmic temperature mean: log(TL/TR)/(1/TR-1/TL).
ratio = deltaT/temperature(2);
if abs(ratio)<1e-7
    meanT = temperature(2)*(1+ratio/2-ratio^2/6+ratio^3/12);
else
    meanT = temperature(1)*log1p(ratio)/ratio;
end
betaJump = deltaT/prod(temperature);
hydraulic = mass*deltaP*mean(1./temperature);
fourier = thermal*deltaT*betaJump;
upwind = mass*p.heatWater_JkgK*(temperature(up)-meanT)*betaJump;
f = p.heatWater_JkgK*(temperature-p.meltingTemperature_K- ...
    temperature.*log(temperature/p.meltingTemperature_K));
mu = f+P;
direct = mass*(mu(1)/temperature(1)-mu(2)/temperature(2))+energy*betaJump;
audit = struct('entropyProduction_WmK',hydraulic+fourier+upwind, ...
    'hydraulic_WmK',hydraulic,'fourier_WmK',fourier,'upwind_WmK',upwind, ...
    'directEntropyProduction_WmK',direct,'massFlux_kgms',mass,'energyFlux_Wm',energy);
end
