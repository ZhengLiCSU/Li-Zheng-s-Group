function [increment,life]=luo2026_fatigue_increment(stressRatio,a,b)
%LUO2026_FATIGUE_INCREMENT Constitutive reproduction of Eqs.25-26.
% Luo et al., Computers and Geotechnics 196 (2026) 108114.
% stressRatio is peak EFFECTIVE tensile stress / original tensile strength.
% It is not a temperature, ice fraction, or tunnel-gasket damage parameter.
% This kernel excludes their finite-element stress redistribution/viscosity.
arguments
    stressRatio double {mustBeReal,mustBeFinite}
    a (1,1) double {mustBeReal,mustBeFinite} = -.051085
    b (1,1) double {mustBeReal,mustBeFinite} = .9929
end
assert(a<0&&all(stressRatio(:)>0&stressRatio(:)<1), ...
    'Luo2026:Domain','Use a<0 and subcritical positive tensile peak ratios.');
logLife=(stressRatio-b)/a;
assert(all(logLife(:)>=0&logLife(:)<300),'Luo2026:SNDomain', ...
    'The supplied S-N fit is outside its finite fatigue-life domain.');
life=10.^logLife;
increment=(1-stressRatio)./life;
end
