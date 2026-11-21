function a = solve_cohesive_ring_mechanics(pg,Ts,m,p,kappaOld,loadingAtHistory,needDerivatives)
%SOLVE_COHESIVE_RING_MECHANICS Exact nonlinear static condensation.
if nargin<6, loadingAtHistory=[]; end
if nargin<7, needDerivatives=true; end
free=m.bodyPressure*pg+m.bodyTemperature*(Ts-m.referenceTemperature_K);
B=m.gapOperator; C=m.gapCompliance; target=B*free; n=numel(target);
gap=(eye(n)+C*m.cohesive.tensileStiffness_Pam)\target;
for iteration=1:40
    local=joint_cohesive_energy(gap,kappaOld,m.cohesive);
    r=gap+C*local.traction_Pa-target;
    tolerance=32*eps*max(1e-12,norm(abs(gap)+abs(C)*abs(local.traction_Pa)+abs(target),inf));
    if norm(r,inf)<=tolerance, break; end
    J=eye(n)+C.*local.tangent_Pam';
    change=-J\r; alpha=1; accepted=false;
    for line=1:30
        trial=gap+alpha*change;
        q=joint_cohesive_energy(trial,kappaOld,m.cohesive);
        rt=trial+C*q.traction_Pa-target;
        if norm(rt,2)<(1-1e-4*alpha)*norm(r,2) || norm(rt,inf)<=tolerance
            gap=trial; accepted=true; break
        end
        alpha=alpha/2;
    end
    assert(accepted,'Cohesive:Mechanics','Nonlinear interface line search failed: residual %.4e, bound %.4e.',norm(r,inf),tolerance);
end
local=joint_cohesive_energy(gap,kappaOld,m.cohesive);
assert(norm(gap+C*local.traction_Pa-target,inf)<=tolerance,'Cohesive:Mechanics','Nonlinear interface equilibrium unresolved.');
u=free-m.bodyTraction*local.traction_Pa;
local=joint_cohesive_energy(B*u,kappaOld,m.cohesive);
% At the maximum-opening kink either tangent belongs to the generalized
% derivative. The coupled solver selects the side consistent with its
% proposed gap direction. Traction, history and energy are never altered.
tie=cohesive_history_ties(B*u,kappaOld,m.cohesive);
if ~isempty(loadingAtHistory)
    soft=-m.cohesive.strength_Pa/(m.cohesive.failureGap_m-m.cohesive.onsetGap_m);
    sec=m.cohesive.strength_Pa*(m.cohesive.failureGap_m-kappaOld)./ ...
        ((m.cohesive.failureGap_m-m.cohesive.onsetGap_m)*max(kappaOld,m.cohesive.onsetGap_m));
    local.tangent_Pam(tie & loadingAtHistory)=soft;
    local.tangent_Pam(tie & ~loadingAtHistory)=sec(tie & ~loadingAtHistory);
    % Return a coherent generalized derivative for all material metadata,
    % including callers that later use fracture-energy sensitivities.
    active=tie & loadingAtHistory;
    local.damageDerivative(tie)=0; local.surfaceDerivative(tie)=0;
    local.damageDerivative(active)=m.cohesive.failureGap_m*m.cohesive.onsetGap_m./ ...
        ((m.cohesive.failureGap_m-m.cohesive.onsetGap_m)*kappaOld(active).^2);
    local.surfaceDerivative(active)=m.cohesive.fractureEnergy_Jm2/(m.cohesive.failureGap_m-m.cohesive.onsetGap_m);
    local.elasticDerivative(tie)=local.traction_Pa(tie)- ...
        .5*m.cohesive.tensileStiffness_Pam*max(B(tie,:)*u,0).^2.*local.damageDerivative(tie);
end
cached=isfield(m,'elasticSolverCache') && ~m.cohesive.enabled && ...
    m.cohesive.tensileStiffness_Pam==m.cohesive.contactStiffness_Pam;
if cached
    K=m.elasticSolverCache.tangent;
else
    K=m.stiffness+B'*diag(m.weight.*local.tangent_Pam)*B;
    [~,flag]=chol(K); assert(flag==0,'Cohesive:Stability','Unsupported unstable equilibrium.');
end
a=struct('u',u,'local',local,'tangent',K);
if needDerivatives
    if cached
        a.pressureDerivative=m.elasticSolverCache.pressureDerivative;
        a.temperatureDerivative=m.elasticSolverCache.temperatureDerivative;
    else
        a.pressureDerivative=K\(B'*diag(m.weight*p.wetAreaFraction));
        a.temperatureDerivative=K\m.thermalLoad;
    end
end
end
