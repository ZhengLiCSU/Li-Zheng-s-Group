function report = cohesive_isothermal_leak_probe(model,p,history,headPressure)
%COHESIVE_ISOTHERMAL_LEAK_PROBE Controlled warm hydraulic reload.
% Hold ALL temperatures and boundary pressures identical across histories.
% This is a diagnostic steady problem, not another transient tunnel year.
arguments
    model (1,1) struct
    p (1,1) struct
    history (:,1) double
    headPressure (1,1) double {mustBePositive} = 2e5
end
n=numel(model.weight); ns=numel(model.capacity); T=model.referenceTemperature_K;
z=[log1p(headPressure/2/p.atmosphericPressure_Pa)*ones(n,1);zeros(2*n,1);T*ones(n+ns,1)];
bc=thermo_interface_boundary(model,T,T,headPressure,false);
ms=model.weight*p.rhoWater*p.aperture0_m;
for iteration=1:40
    [r,J]=cohesive_threephase_residual(z,z,1,model,bc,p,history);
    if norm(r(1:n),inf)<1e-12, break; end
    [~,d]=cohesive_threephase_state(z,model,p,history);
    % The old state equals the trial state here: subtract the transient
    % storage derivative to obtain the steady hydraulic Jacobian.
    A=J(1:n,1:n)-(model.weight./ms).*d.mass(:,1:n);
    change=-A\r(1:n); alpha=1; accepted=false;
    for line=1:30
        trial=z; trial(1:n)=z(1:n)+alpha*change;
        rt=cohesive_threephase_residual(trial,trial,1,model,bc,p,history);
        if norm(rt(1:n),2)<(1-1e-4*alpha)*norm(r(1:n),2)
            z=trial; accepted=true; break;
        end
        alpha=alpha/2;
    end
    assert(accepted,'Cohesive:Probe','Steady hydraulic reload did not converge.');
end
[r,~,a]=cohesive_threephase_residual(z,z,1,model,bc,p,history);
assert(norm(r(1:n),inf)<1e-12,'Cohesive:Probe','Steady mass balance failed.');
assert(all(a.state.cohesive.kappa_m<=max(history,model.cohesive.onsetGap_m)+1e-12), ...
    'Cohesive:Probe','Reload caused additional fracture; not a fixed-damage comparison.');
outlets=model.faces(:,2)==0;
report=struct('rate_Lmday',sum(a.faceFlux(outlets,1))*86400*1000/p.rhoWater, ...
    'jointRate_Lmday',a.faceFlux(outlets,1)*86400*1000/p.rhoWater, ...
    'pressure_Pa',a.state.pressure,'gap_m',a.state.x(:,1), ...
    'temperature_K',T,'reservoirPressure_Pa',headPressure, ...
    'maxDamage',max(a.state.cohesive.damage),'massResidual_kgms',ms.*r(1:n), ...
    'state',z,'history',history,'scope','Steady fully warm reload with fixed previous fracture history');
end
