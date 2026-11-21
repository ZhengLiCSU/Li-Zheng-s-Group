function c = joint_cohesive_parameters(overrides)
%JOINT_COHESIVE_PARAMETERS Illustrative mode-I bonded interface, NOT EPDM.
% Strength and fracture energy require interface-specific calibration.
% Entire irreversible fracture work is stored as surface energy in this
% ideal energetic limit; no ad hoc damage-to-heat conversion is introduced.
arguments
    overrides (1,1) struct = struct()
end
c=struct('tensileStiffness_Pam',1e12,'contactStiffness_Pam',1e12, ...
    'strength_Pa',2e6,'fractureEnergy_Jm2',200,'enabled',true);
names=fieldnames(overrides);
for k=1:numel(names)
    assert(isfield(c,names{k}),'Cohesive:Parameter','Unknown material field.');
    c.(names{k})=overrides.(names{k});
end
values=[c.tensileStiffness_Pam,c.contactStiffness_Pam,c.strength_Pa,c.fractureEnergy_Jm2];
assert(all(isfinite(values)&values>0),'Cohesive:Parameter','Positive finite material values required.');
c.onsetGap_m=c.strength_Pa/c.tensileStiffness_Pam;
c.failureGap_m=2*c.fractureEnergy_Jm2/c.strength_Pa;
assert(c.failureGap_m>c.onsetGap_m,'Cohesive:Parameter','Fracture work must exceed the elastic onset energy.');
assert(islogical(c.enabled)&&isscalar(c.enabled),'Cohesive:Parameter','enabled must be scalar logical.');
end
