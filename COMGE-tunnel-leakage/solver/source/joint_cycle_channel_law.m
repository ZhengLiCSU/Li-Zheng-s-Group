function [law,info]=joint_cycle_channel_law(base,severity,c)
%JOINT_CYCLE_CHANNEL_LAW Existing-channel coarsening at CONSTANT void volume.
% Uniform radius dilation lambda and count-density reduction lambda^-2 give
% Cref'=lambda^2*Cref, rQ'=lambda*rQ, V'=V. Fixed lognormal shape/path.
% This is an explicit research hypothesis, not calibrated gasket fracture.
c=joint_cycle_damage_defaults(c);law=base;nj=numel(base.jointAnglesRad);
% Accepted-history blocks always contain one conductance per joint.
if nj>1&&isscalar(law.jointConductance)
    law.jointConductance=repmat(law.jointConductance,nj,1);
end
D=expand(severity,nj);alpha=expand(c.radiusGrowthAtUnitDamage,nj);
assert(all(D>=0&D<=1),'JointCycle:Severity','Severity must lie in [0,1].');
lambda=1+alpha.*D;
assert(all(isfinite(lambda.^2)),'JointCycle:Geometry','Coarsening geometry overflows.');
if any(lambda~=1)
    law.jointConductance=expand(base.jointConductance,nj).*lambda.^2;
    law.jointFreeze.flowWeightedRadius_m=expand(base.jointFreeze.flowWeightedRadius_m,nj).*lambda;
end
info=struct('severity',D,'radiusMultiplier',lambda, ...
    'channelDensityMultiplier',lambda.^(-2),'conductanceMultiplier',lambda.^2, ...
    'connectedVolumeMultiplier',ones(nj,1),'evolution','accepted_local_freeze_thaw_coarsening', ...
    'scope','Illustrative existing-channel redistribution; constant connected volume; no seal breakthrough or mechanical softening.');
end
function v=expand(v,n)
assert(isnumeric(v)&&isreal(v)&&isvector(v)&&all(isfinite(v(:))) ...
    &&(isscalar(v)||numel(v)==n),'JointCycle:Shape','Require finite scalar or joint vector.');
v=double(v(:));if isscalar(v),v=repmat(v,n,1);end
end
