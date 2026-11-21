function c=joint_cycle_damage_defaults(user)
%JOINT_CYCLE_DAMAGE_DEFAULTS Optional, uncalibrated hydraulic-coarsening law.
% The cyclic rock S-N parameters of Luo2026 are deliberately NOT used here.
c=struct('enabled',false,'initialSeverity',0,'doseRate',.05, ...
    'doseExponent',1,'radiusGrowthAtUnitDamage',.5, ...
    'freezeFraction',.2,'thawFraction',.02,'minimumFrozenDays',0);
if nargin>0
    assert(isstruct(user)&&isscalar(user),'JointCycle:Settings','Settings must be scalar struct.');
    for key=fieldnames(user)'
        assert(isfield(c,key{1}),'JointCycle:Settings','Unknown setting: %s.',key{1});
        c.(key{1})=user.(key{1});
    end
end
assert((isnumeric(c.enabled)||islogical(c.enabled))&&isscalar(c.enabled) ...
    &&isreal(c.enabled)&&isfinite(c.enabled)&&ismember(c.enabled,[0 1]), ...
    'JointCycle:Settings','enabled must be scalar 0/1.');
for key={'initialSeverity','doseRate','radiusGrowthAtUnitDamage'}
    v=c.(key{1});
    assert(isnumeric(v)&&isreal(v)&&isvector(v)&&all(isfinite(v(:)))&&all(v(:)>=0), ...
        'JointCycle:Settings','%s must be finite nonnegative scalar or joint vector.',key{1});
end
assert(all(c.initialSeverity(:)<=1),'JointCycle:Settings','Initial severity must be in [0,1].');
for key={'doseExponent','freezeFraction','thawFraction','minimumFrozenDays'}
    v=c.(key{1});
    assert(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v), ...
        'JointCycle:Settings','%s must be a finite scalar.',key{1});
end
assert(c.doseExponent>0&&c.thawFraction>=0&&c.freezeFraction<=1 ...
    &&c.freezeFraction>c.thawFraction&&c.minimumFrozenDays>=0, ...
    'JointCycle:Settings','Require 0<=thaw<freeze<=1, positive exponent and nonnegative duration.');
end
