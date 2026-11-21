function report=audit_joint_cycle_damage(out)
%AUDIT_JOINT_CYCLE_DAMAGE Replay saved accepted ice-state cycles and geometry.
% Independent exponential moment/volume identities supplement memory replay.
% This audits numerical records, not the empirical validity of coarsening.
c=out.jointCycleDamage;assert(c.enabled,'JointCycleAudit:Disabled','No evolving case.');
p=out.parameters;h=c.history;n=size(out.history,1);base=joint_microchannel_damage(p.leakage);
nj=numel(base.jointAnglesRad);nq=size(c.initialMemory.lastIce,2);
assert(isequal(size(h),[n 49])&&isequal(size(c.sampleIceHistory),[n nj*nq]) ...
    &&all(isfinite(h(:)))&&all(isfinite(c.sampleIceHistory(:))), ...
    'JointCycleAudit:History','Incomplete/nonfinite cycle records.');
alpha=expand(c.settings.radiusGrowthAtUnitDamage,nj);
rate=expand(c.settings.doseRate,nj);D0=expand(c.settings.initialSeverity,nj);
memory=c.initialMemory;maxReplay=0;maxConductance=0;maxVolume=0;
[~,~,inventory]=joint_channel_inventory(zeros(nj,1),base.jointFreeze, ...
    base.jointStorage,base.jointConductance,p.liningThickness,p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
c0=expand(base.jointConductance,nj);active=c0>0;
for k=1:n
    maxReplay=max(maxReplay,max(abs(h(k,2:9)'-memory.targetSeverity)));
    ice=reshape(c.sampleIceHistory(k,:),nq,nj)';
    memory=joint_freeze_thaw_memory(memory,ice,h(k,1),active,c.settings);
    maxReplay=max(maxReplay,max(abs(h(k,10:17)'-memory.targetSeverity)));
    Dcheck=D0+(1-D0).*(-expm1(-rate.*memory.cumulativeDose));
    maxReplay=max(maxReplay,max(abs(mean(Dcheck,2)-h(k,10:17)')));
    maxReplay=max(maxReplay,max(abs(mean(memory.cycleCount,2)-h(k,18:25)')));
    lambda=1+alpha.*h(k,2:9)';expected=c0.*lambda.^2;
    maxConductance=max(maxConductance,max(abs(h(k,42:49)'-expected)));
    maxVolume=max(maxVolume,max(abs(h(k,34:41)'-inventory.waterVolume_m2)));
end
report=struct('passed',maxReplay<1e-12&&maxConductance<1e-20&&maxVolume<1e-15 ...
    &&isequaln(memory,c.finalMemory)&&isequal(h(:,1),out.history(:,1)), ...
    'maxReplayDifference',maxReplay,'maxConductanceDifference_m_s',maxConductance, ...
    'maxVolumeDifference_m2',maxVolume,'finalMemoryExact',isequaln(memory,c.finalMemory), ...
    'acceptedSteps',n,'scope','Accepted-state memory and constitutive inventory audit; no experimental validation.');
end
function v=expand(v,n)
v=v(:);if isscalar(v),v=repmat(v,n,1);end
end
