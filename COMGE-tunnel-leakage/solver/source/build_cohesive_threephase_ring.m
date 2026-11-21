function m = build_cohesive_threephase_ring(S,nq,p,c,options)
%BUILD_COHESIVE_THREEPHASE_RING Nonlinear interface, unchanged Q4 bulk.
arguments
    S (1,1) struct = make_thm_verification_mesh(16)
    nq (1,1) double {mustBePositive,mustBeInteger} = 3
    p (1,1) struct = joint_threephase_parameters()
    c (1,1) struct = joint_cohesive_parameters()
    options (1,1) struct = struct()
end
options.normalRestraint_Pam=c.tensileStiffness_Pam;
m=build_thermo_interface_ring(S,nq,p,options); B=m.gapOperator; w=m.weight;
m.stiffness=m.stiffness-B'*diag(w*c.tensileStiffness_Pam)*B;
m.stiffness=(m.stiffness+m.stiffness')/2;
% Condense the linear bulk only; the nonlinear interface is solved each
% constitutive trial with the PREVIOUS ACCEPTED history held fixed.
m.bodyPressure=m.stiffness\(B'*diag(w*p.wetAreaFraction));
m.bodyTemperature=m.stiffness\m.thermalLoad;
m.bodyTraction=m.stiffness\(B'*diag(w));
m.gapCompliance=B*m.bodyTraction; m.cohesive=c;
m.scope='Q4 bulk + energetic mode-I cohesive contact + pure-water three-phase slit; no soil or calibrated gasket.';
% Global convexity over ALL pieces guarantees a unique supported equilibrium.
soft=0;
if c.enabled, soft=-c.strength_Pa/(c.failureGap_m-c.onsetGap_m); end
lower=m.stiffness+B'*diag(w*min(soft,0))*B;
[~,flag]=chol(lower);
assert(flag==0,'Cohesive:Stability', ...
    'Supported tangent can lose positive definiteness. Arc-length/dynamics needed; do not clip the softening tangent.');
end
