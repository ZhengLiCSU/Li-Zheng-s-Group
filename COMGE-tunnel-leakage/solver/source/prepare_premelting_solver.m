function m = prepare_premelting_solver(m,p)
%PREPARE_PREMELTING_SOLVER Exact reusable tangent for a linear elastic joint.
% Cache only when damage is disabled AND tensile/contact stiffnesses match.
% Geometry, stiffness, thermal force, weights and wet fraction are keys.
% Contact with unequal stiffness or any damage retains the original solve.
c=m.cohesive;
eligible=~c.enabled && c.tensileStiffness_Pam==c.contactStiffness_Pam;
if ~eligible
    if isfield(m,'elasticSolverCache'), m=rmfield(m,'elasticSolverCache'); end
    return
end
key=struct('stiffness',m.stiffness,'gapOperator',m.gapOperator, ...
    'weight',m.weight,'thermalLoad',m.thermalLoad,'wetAreaFraction',p.wetAreaFraction, ...
    'tensileStiffness_Pam',c.tensileStiffness_Pam,'contactStiffness_Pam',c.contactStiffness_Pam);
if isfield(m,'elasticSolverCache') && isequal(m.elasticSolverCache.key,key), return; end
B=m.gapOperator;
K=m.stiffness+B'*diag(m.weight*c.tensileStiffness_Pam)*B;
[~,flag]=chol(K); assert(flag==0,'Cohesive:Stability','Unsupported unstable elastic equilibrium.');
cache=struct('key',key,'tangent',K, ...
    'pressureDerivative',K\(B'*diag(m.weight*p.wetAreaFraction)), ...
    'temperatureDerivative',K\m.thermalLoad);
m.elasticSolverCache=cache;
end
