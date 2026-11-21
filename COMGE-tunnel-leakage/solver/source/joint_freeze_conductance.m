function [factor,dFactor_dT,diagnostic] = joint_freeze_conductance(Tprofile,settings)
%JOINT_FREEZE_CONDUCTANCE Reversible ice/viscosity conductance multiplier.
%   TPROFILE is [number of joints x number of radial midpoint samples] in C.
%   FACTOR multiplies conductance calibrated at referenceTemperature_C.
%   DFACTORT is the derivative of each joint's factor with respect to each
%   temperature sample, with the same shape as TPROFILE.
%
%   This is a locally homogenized, interconnected microchannel-network
%   approximation. Local mobilities are combined as hydraulic resistances
%   in series with equal midpoint length weights. It is not the exact result
%   for isolated parallel capillaries of persistent radius through the whole
%   lining. Local thermal equilibrium with the lining is assumed; channel
%   latent storage, concentrated advective warming and permanent damage are
%   not represented by this reduced constitutive law.
%
%   For cylindrical capillaries, Gibbs-Thomson gives
%       rcrit = 2*sigma_iw*Tm/(rho_i*Lf*(Tbulk-T)).
%   Larger pores freeze first. If ln(r) ~ N(m,s^2), the r^4 flow-weighted
%   median is rQ = exp(m+4*s^2). The remaining capillary conductance is
%       Fice = E[r^4*I(r<rcrit)]/E[r^4]
%            = Phi(log(rcrit/rQ)/s),       T < Tbulk,
%            = 1,                        T >= Tbulk.
%   This truncated lognormal moment is a reduced-model derivation, not an
%   experimentally calibrated tunnel-gasket freezing curve. Small surviving
%   channels supply a smooth tail; no arbitrary residual-flow floor is used.
%
%   With midpoint weights 1/nq and muRef=mu(referenceTemperature_C),
%       Rk = mu(Tk)/(muRef*Fice(Tk)),
%       factor = 1/mean(Rk),
%       dfactor/dTk = -factor * (Rk/sum(R)) * dlog(Rk)/dTk.
%   Evaluate Fice and the series resistance in log space to avoid premature
%   loss of tiny conductance tails. Genuine floating-point underflow to zero
%   is possible at extreme blocking; it is not a prescribed zero-flow rule.
%
%   Sources:
%   Luo et al. (2023), JGR: Solid Earth, doi:10.1029/2022JB025254,
%       https://doi.org/10.1029/2022JB025254 (Gibbs-Thomson/geometry/constants).
%   Wang and Vanapalli (2025), WRR, doi:10.1029/2024WR039453,
%       https://doi.org/10.1029/2024WR039453 (capillary and film flow limits).
%   IAPWS SR6-08(2011), Eq. 7 and Table 5, liquid-water viscosity at 0.1 MPa,
%       https://www.iapws.org/relguide/LiquidWater.pdf
%       Valid -20 to 110 C, including metastable liquid water below 0 C.
%   No adsorption-film contribution is calibrated here; this law should not
%   be interpreted as predicting all film leakage in a frozen gasket.
if nargin < 2 || isempty(settings)
    settings = joint_freeze_defaults();
else
    settings = joint_freeze_defaults(settings);
end
assert(isnumeric(Tprofile) && isreal(Tprofile) && ismatrix(Tprofile) && ...
    all(isfinite(Tprofile(:))), ...
    'joint_freeze_conductance:InvalidTemperature', ...
    'Tprofile must be a finite real numeric joint-by-sample matrix.');
Tprofile = double(Tprofile);
[nj,nq] = size(Tprofile);
[radiusProfile,geometryMobility]=joint_radius_profile(settings,nj,nq);
assert(nq > 0 || nj == 0, ...
    'joint_freeze_conductance:EmptyPath', ...
    'Each joint needs at least one radial temperature sample.');
logicalNames = {'enabled','iceEnabled','viscosityEnabled','heatExchangeEnabled'};
for k = 1:numel(logicalNames)
    value = settings.(logicalNames{k});
    assert(isscalar(value) && (islogical(value) || isnumeric(value)) && ...
        isreal(value) && isfinite(value) && any(value == [0 1]), ...
        'joint_freeze_conductance:InvalidSwitch', ...
        '%s must be scalar logical or numeric 0/1.',logicalNames{k});
end
assert(isnumeric(settings.bulkFreezingTemperature_C) && ...
    isreal(settings.bulkFreezingTemperature_C) && ...
    isscalar(settings.bulkFreezingTemperature_C) && ...
    isfinite(settings.bulkFreezingTemperature_C), ...
    'joint_freeze_conductance:InvalidBulkFreezingTemperature', ...
    'bulkFreezingTemperature_C must be a finite real scalar.');
assert(isnumeric(settings.quadraturePoints) && ...
    isreal(settings.quadraturePoints) && isscalar(settings.quadraturePoints) && ...
    isfinite(settings.quadraturePoints) && settings.quadraturePoints >= 1 && ...
    settings.quadraturePoints == fix(settings.quadraturePoints), ...
    'joint_freeze_conductance:InvalidQuadraturePoints', ...
    'quadraturePoints must be a positive integer.');

factor = ones(nj,1);
dFactor_dT = zeros(size(Tprofile));
iceFraction = ones(size(Tprofile));
mobility = ones(size(Tprofile));
frozenLength = mean(Tprofile < settings.bulkFreezingTemperature_C,2);
% This diagnostic is the geometric below-bulk-freezing path fraction. It
% does not assert that all water in that part of the path has become ice.
if nj == 0
    minTemperature = zeros(0,1);
    maxTemperature = zeros(0,1);
else
    minTemperature = min(Tprofile,[],2);
    maxTemperature = max(Tprofile,[],2);
end
if ~settings.enabled || nj == 0
    diagnostic = makeDiagnostic(minTemperature,maxTemperature, ...
        iceFraction,mobility,factor,frozenLength);
    return
end

logIce = zeros(size(Tprofile));
dLogIce_dT = zeros(size(Tprofile));
if settings.iceEnabled
    radiusQ = jointParameter(settings.flowWeightedRadius_m,nj, ...
        'flowWeightedRadius_m');
    logStd = jointParameter(settings.logRadiusStd,nj,'logRadiusStd');
    positiveNames = {'interfacialEnergy_J_m2','iceDensity_kg_m3', ...
        'latentHeat_J_kg','meltingTemperature_K'};
    for k = 1:numel(positiveNames)
        value = settings.(positiveNames{k});
        assert(isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value) && value > 0, ...
            'joint_freeze_conductance:InvalidPhysicalConstant', ...
            '%s must be a positive finite real scalar.',positiveNames{k});
    end
    % Work with log(beta) so exceptionally small/large constants do not
    % overflow the dimensional critical-radius product or quotient.
    logBeta = log(2) + log(settings.interfacialEnergy_J_m2) + ...
        log(settings.meltingTemperature_K) - log(settings.iceDensity_kg_m3) - ...
        log(settings.latentHeat_J_kg);
    depression = settings.bulkFreezingTemperature_C - Tprofile;
    cold = depression > 0;
    if any(cold(:))
        logRadiusQ = log(radiusProfile);
        stdGrid = repmat(logStd,1,nq);
        z = (logBeta - log(depression(cold)) - logRadiusQ(cold)) ./ stdGrid(cold);
        [logPhi,mills] = normalLogCDF(z);
        logIce(cold) = logPhi;
        dLogIce_dT(cold) = mills ./ stdGrid(cold) ./ depression(cold);
    end
    iceFraction = exp(logIce);
end

logViscosityRatio = zeros(size(Tprofile));
dLogMu_dT = zeros(size(Tprofile));
if settings.viscosityEnabled
    assert(isnumeric(settings.referenceTemperature_C) && ...
        isreal(settings.referenceTemperature_C) && ...
        isscalar(settings.referenceTemperature_C) && ...
        isfinite(settings.referenceTemperature_C), ...
        'joint_freeze_conductance:InvalidReferenceTemperature', ...
        'referenceTemperature_C must be a finite real scalar.');
    [mu,dMu_dT] = waterViscosity(Tprofile);
    muReference = waterViscosity(settings.referenceTemperature_C);
    logViscosityRatio = log(mu) - log(muReference);
    dLogMu_dT = dMu_dT ./ mu;
end
logResistance = logViscosityRatio - logIce - log(geometryMobility);
mobility = exp(-logResistance);

maxLogResistance = max(logResistance,[],2);
finiteRows = isfinite(maxLogResistance);
factor(~finiteRows) = 0;
if any(finiteRows)
    resistanceScaled = exp(logResistance(finiteRows,:) - ...
        maxLogResistance(finiteRows));
    sumScaled = sum(resistanceScaled,2);
    logMeanResistance = maxLogResistance(finiteRows) + log(sumScaled) - log(nq);
    factor(finiteRows) = exp(-logMeanResistance);
    derivativeRows = find(finiteRows & factor > 0);
    if ~isempty(derivativeRows)
        weights = exp(logResistance(derivativeRows,:) - ...
            maxLogResistance(derivativeRows));
        weights = weights ./ sum(weights,2);
        dLogResistance_dT = dLogMu_dT(derivativeRows,:) - ...
            dLogIce_dT(derivativeRows,:);
        % A vanishing local weight contributes zero even if an extreme
        % parameter choice makes its log derivative overflow.
        contributions = weights .* dLogResistance_dT;
        contributions(weights == 0) = 0;
        dFactor_dT(derivativeRows,:) = -factor(derivativeRows) .* contributions;
    end
end
diagnostic = makeDiagnostic(minTemperature,maxTemperature, ...
    iceFraction,mobility,factor,frozenLength);
end

function value = jointParameter(value,nj,name)
assert(isnumeric(value) && isreal(value) && isvector(value) && ...
    all(isfinite(value(:))) && all(value(:) > 0) && ...
    (isscalar(value) || numel(value) == nj), ...
    'joint_freeze_conductance:InvalidJointParameter', ...
    '%s must be positive finite scalar or one value per joint.',name);
value = double(value(:));
if isscalar(value)
    value = repmat(value,nj,1);
end
end

function [logPhi,mills] = normalLogCDF(z)
% log(Phi(z)) and phi(z)/Phi(z), stable on the negative normal tail.
logPhi = zeros(size(z));
mills = zeros(size(z));
negative = z < 0;
if any(negative)
    scaledTail = erfcx(-z(negative)/sqrt(2));
    logPhi(negative) = log(0.5) + log(scaledTail) - 0.5*z(negative).^2;
    mills(negative) = sqrt(2/pi) ./ scaledTail;
end
nonnegative = ~negative;
if any(nonnegative)
    logPhi(nonnegative) = log1p(-0.5*erfc(z(nonnegative)/sqrt(2)));
    mills(nonnegative) = exp(-0.5*z(nonnegative).^2 - ...
        0.5*log(2*pi) - logPhi(nonnegative));
end
end

function [mu,dMu_dT] = waterViscosity(temperature_C)
% Official four-term ambient-pressure IAPWS equation and its derivative.
assert(isnumeric(temperature_C) && isreal(temperature_C) && ...
    all(isfinite(temperature_C(:))) && ...
    all(temperature_C(:) >= -20 & temperature_C(:) <= 110), ...
    'joint_freeze_conductance:ViscosityTemperatureRange', ...
    'The IAPWS liquid-water viscosity law requires -20 <= T <= 110 C.');
tau = (double(temperature_C) + 273.15)/300;
coefficients = [280.68,511.45,61.131,0.45903];
exponents = [-1.9,-7.7,-19.6,-40];
mu = zeros(size(tau));
dMu_dT = zeros(size(tau));
for k = 1:4
    mu = mu + 1e-6*coefficients(k)*tau.^exponents(k);
    dMu_dT = dMu_dT + (1e-6*coefficients(k)*exponents(k)/300) * ...
        tau.^(exponents(k)-1);
end
end

function diagnostic = makeDiagnostic(tmin,tmax,iceFraction,mobility,factor,frozenLength)
diagnostic = struct('minTemperature_C',tmin, ...
    'maxTemperature_C',tmax, ...
    'iceConductanceFraction',iceFraction, ...
    'mobility',mobility, ...
    'factor',factor, ...
    'effectiveFrozenLengthFraction',frozenLength);
end
