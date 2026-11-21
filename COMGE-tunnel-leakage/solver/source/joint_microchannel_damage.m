function [law,info] = joint_microchannel_damage(leakage)
%JOINT_MICROCHANNEL_DAMAGE Resolve a fixed residual joint-channel geometry.
% Damage is a normalized HYDRAULIC geometry coordinate, not 1-E/E0 and not
% a number of freeze-thaw cycles. It describes existing connected interface
% channels, not concrete-body cracks or molecular pores inside EPDM rubber.
% For a fixed lognormal shape and path length, dilating every radius by
% lambda=1+radiusGrowthAtUnitDamage*D and channel density by eta gives
%     rQ(D)=lambda*rQ0; Cref(D)=eta*lambda^4*Cref0;
%     connectedVolume(D)=eta*lambda^2*connectedVolume0.
% Both conductance and the r^2-weighted water/ice inventory MUST use LAW.
% The changed radius also shifts Gibbs-Thomson freezing, so conductance at
% subzero temperature need not increase monotonically with damage.
% Zero baseline connectivity stays zero. No new breakthrough is invented.
% Default disabled; no uncalibrated cycle-to-damage evolution is applied.
arguments
    leakage (1,1) struct
end
settings=struct('enabled',false,'severity',0, ...
    'radiusGrowthAtUnitDamage',0.5,'channelDensityGrowthAtUnitDamage',0);
if isfield(leakage,'jointDamage')
    supplied=leakage.jointDamage;
    assert(isstruct(supplied)&&isscalar(supplied),'JointDamage:Settings', ...
        'jointDamage must be a scalar structure.');
    keys=fieldnames(supplied);
    for k=1:numel(keys)
        assert(isfield(settings,keys{k}),'JointDamage:Settings','Unknown damage setting: %s.',keys{k});
        settings.(keys{k})=supplied.(keys{k});
    end
end
assert((isnumeric(settings.enabled)||islogical(settings.enabled)) ...
    &&isscalar(settings.enabled)&&isreal(settings.enabled)&&isfinite(settings.enabled) ...
    &&ismember(settings.enabled,[0 1]),'JointDamage:Switch','enabled must be logical 0/1.');
law=leakage;
nj=numel(leakage.jointAnglesRad);
info=struct('enabled',logical(settings.enabled),'settings',settings, ...
    'severity',zeros(nj,1),'radiusMultiplier',ones(nj,1), ...
    'channelDensityMultiplier',ones(nj,1),'referenceConductanceMultiplier',ones(nj,1), ...
    'connectedVolumeMultiplier',ones(nj,1),'evolution','fixed_residual_geometry', ...
    'scope','Existing joint interface channels; illustrative geometry mapping, no gasket/contact failure or cumulative damage calibration.');
if ~settings.enabled
    return
end
assert(nj>0&&leakage.jointFreeze.enabled&&leakage.jointStorage.enabled, ...
    'JointDamage:Inventory','Damage requires joint thermal paths and enabled liquid/ice inventory.');
D=jointVector(settings.severity,nj,'Severity');
alpha=jointVector(settings.radiusGrowthAtUnitDamage,nj,'RadiusGrowth');
beta=jointVector(settings.channelDensityGrowthAtUnitDamage,nj,'DensityGrowth');
assert(all(D>=0&D<=1)&&all(alpha>=0)&&all(beta>=0), ...
    'JointDamage:Range','severity must be in [0,1]; geometry-growth parameters must be nonnegative.');
c0=jointVector(leakage.jointConductance,nj,'Conductance');
rq0=jointVector(leakage.jointFreeze.flowWeightedRadius_m,nj,'Radius');
assert(all(c0>=0)&&all(rq0>0),'JointDamage:BaseGeometry','Supply nonnegative Cref and positive rQ.');
lambda=1+alpha.*D;
eta=1+beta.*D;
cScale=eta.*lambda.^4;
vScale=eta.*lambda.^2;
assert(all(isfinite(cScale))&&all(isfinite(vScale)) ...
    &&all(isfinite(c0.*cScale))&&all(isfinite(rq0.*lambda)), ...
    'JointDamage:Overflow','Residual channel geometry overflows; revise physical inputs.');
% Avoid altering scalar/vector shapes when the physical change is exactly zero.
if any(cScale~=1)
    law.jointConductance=c0.*cScale;
end
if any(lambda~=1)
    law.jointFreeze.flowWeightedRadius_m=rq0.*lambda;
end
info.severity=D;
info.radiusMultiplier=lambda;
info.channelDensityMultiplier=eta;
info.referenceConductanceMultiplier=cScale;
info.connectedVolumeMultiplier=vScale;
info.referenceConductance_m_s=c0.*cScale;
info.flowWeightedRadius_m=rq0.*lambda;
end

function value=jointVector(value,nj,name)
assert(isnumeric(value)&&isreal(value)&&isvector(value)&&all(isfinite(value(:))) ...
    &&(isscalar(value)||numel(value)==nj),['JointDamage:' name], ...
    '%s must be finite scalar or one value per joint.',name);
value=double(value(:));
if isscalar(value),value=repmat(value,nj,1);end
end
