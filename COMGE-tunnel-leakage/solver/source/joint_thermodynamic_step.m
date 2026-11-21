function [next,info] = joint_thermodynamic_step(old,temperature,dt,settings,p)
%JOINT_THERMODYNAMIC_STEP Implicit conservative interface step with audits.
% Exact analytic tangent, primal active set for the phase bounds, and a
% feasible line search. Pressure is never clipped. A failed physical domain
% or nonlinear solve is returned as an error, not disguised as zero leakage.
arguments
    old (:,4) double {mustBeReal,mustBeFinite}
    temperature (1,1) double {mustBePositive,mustBeFinite}
    dt (1,1) double {mustBePositive,mustBeFinite}
    settings (1,1) struct
    p (1,1) struct = joint_thermodynamic_parameters()
end
settings = defaults(settings);
oldState = joint_helmholtz_state(old,p);
n = size(old,1); dx = settings.length_m/n;
scaling = repmat([p.aperture0_m/p.wetAreaFraction;1e6;1],n,1);
z = [old(:,1),oldState.pressure_Pa,old(:,3)./old(:,2)];
if p.phaseViscosity_Jm2skg2==0
    if temperature>=p.meltingTemperature_K
        % A liquid predictor reconstructs mass in the iterate only. No
        % physical history is reset; the complete mass balance is solved.
        preload = p.normalStiffness_Pam*p.referenceGap_m+p.externalTraction_Pa;
        predictedGap = (p.wetAreaFraction*z(:,2)+preload)/p.normalStiffness_Pam;
        contact = predictedGap<0;
        predictedGap(contact) = (p.wetAreaFraction*z(contact,2)+preload)/ ...
            (p.normalStiffness_Pam+p.contactStiffness_Pam);
        z(:,1) = predictedGap; z(:,3) = 0;
    else
        try
            [predicted,predictedState] = joint_closed_equilibrium(old(:,2),temperature,p);
            z = [predicted(:,1),predictedState.pressure_Pa,predicted(:,3)./predicted(:,2)];
        catch exception
            if ~startsWith(exception.identifier,'JointThermo:'), rethrow(exception); end
            % A closed-cell predictor can leave the saturated domain on
            % thawing after drainage. Keep the valid old iterate in that case.
        end
    end
end
iterations = 0;
for iteration = 1:settings.maxIterations
    [residual,tangent] = joint_interface_primary_residual(z,old,temperature,dt,settings,p);
    normResidual = norm(residual,inf);
    if settings.trace && mod(iteration,10)==1
        fprintf('Joint Newton %d: R %.5g, phase %.5g..%.5g\n', ...
            iteration,normResidual,min(z(:,3)),max(z(:,3)));
    end
    if normResidual<=settings.tolerance
        iterations = iteration-1;
        break
    end
    update = newtonDirection(tangent,residual,z,scaling);
    update = reshape(scaling.*update,3,n)';
    alpha = feasibleStep(z,update,p);
    accepted = false;
    for backtrack = 1:settings.maxBacktracks
        trial = z+alpha*update;
        trial(:,3) = roundoffBounds(trial(:,3));
        if alpha>0 && admissible(trial,temperature,p)
            trialResidual = joint_interface_primary_residual(trial,old,temperature,dt,settings,p);
            if norm(trialResidual,2)<=(1-1e-4*alpha)*norm(residual,2) || ...
                    norm(trialResidual,inf)<=settings.tolerance
                z = trial; accepted = true; break
            end
        end
        alpha = alpha/2;
    end
    assert(accepted,'JointInterface:LineSearch', ...
        'No feasible step (iteration %d, residual %.3g, alpha %.3g, fraction %.4g..%.4g).', ...
        iteration,normResidual,alpha,min(z(:,3)),max(z(:,3)));
    iterations = iteration;
end
[residual,~,flow,~,x] = joint_interface_primary_residual(z,old,temperature,dt,settings,p);
assert(norm(residual,inf)<=settings.tolerance,'JointInterface:NonConvergence', ...
    'Interface Newton residual %.3g exceeds %.3g.',norm(residual,inf),settings.tolerance);
next = [x,repmat(temperature,n,1)];
s = joint_helmholtz_state(next,p,true,z(:,2));
oldAtNewT = old; oldAtNewT(:,4) = temperature;
oldThermostated = joint_helmholtz_state(oldAtNewT,p);
increment = next(:,1:3)-old(:,1:3);
massDefect = dx*sum(increment(:,2))-dt*(flow.massFlux_kgms(1)-flow.massFlux_kgms(end));
cellMassDefect = dx*increment(:,2)+dt*diff(flow.massFlux_kgms);
chemicalWork = dt*(flow.reservoirPotential_Jkg(1)*flow.massFlux_kgms(1)- ...
    flow.reservoirPotential_Jkg(2)*flow.massFlux_kgms(end));
darcyDissipation = dt*sum(flow.massFlux_kgms.*flow.potentialDifference_Jkg);
phaseDissipation = -dx*sum(s.phaseAffinity_Jkg.*increment(:,3));
deltaFreeEnergy = dx*sum(s.freeEnergy_Jm2-oldThermostated.freeEnergy_Jm2);
convexRemainder = dx*sum(sum(s.gradient(:,1:3).*increment,2)- ...
    (s.freeEnergy_Jm2-oldThermostated.freeEnergy_Jm2));
energyIdentityDefect = deltaFreeEnergy-chemicalWork+darcyDissipation+ ...
    phaseDissipation+convexRemainder;
sl = p.heatWater_JkgK*log(temperature/p.meltingTemperature_K);
hBoundary = flow.reservoirPotential_Jkg+temperature*sl;
advectiveEnergy = dt*(hBoundary(1)*flow.massFlux_kgms(1)-hBoundary(2)*flow.massFlux_kgms(end));
mechanicalWork = p.externalTraction_Pa*dx*sum(increment(:,1));
thermostatHeat = dx*sum(s.internalEnergy_Jm2-oldThermostated.internalEnergy_Jm2)- ...
    advectiveEnergy-mechanicalWork;
entropyProduced = dx*sum(s.entropy_Jm2K-oldThermostated.entropy_Jm2K)- ...
    thermostatHeat/temperature-dt*sl*(flow.massFlux_kgms(1)-flow.massFlux_kgms(end));
temperatureJumpHeat = dx*sum(oldThermostated.internalEnergy_Jm2-oldState.internalEnergy_Jm2);
temperatureJumpEntropy = dx*sum(oldThermostated.entropy_Jm2K-oldState.entropy_Jm2K);
energyScale = max([1,abs(chemicalWork),abs(darcyDissipation),abs(deltaFreeEnergy)]);
assert(max(abs(cellMassDefect))<=settings.massAuditTolerance_kgm, ...
    'JointInterface:MassAudit','Local mass defect exceeds the audit threshold.');
assert(abs(energyIdentityDefect)<=settings.energyAuditRelativeTolerance*energyScale, ...
    'JointInterface:EnergyAudit','Free-energy identity defect exceeds the audit threshold.');
assert(min([darcyDissipation,phaseDissipation,convexRemainder,temperature*entropyProduced])>= ...
    -settings.energyAuditRelativeTolerance*energyScale, ...
    'JointInterface:DissipationAudit','A negative dissipation or convexity defect was detected.');
info = struct('iterations',iterations,'scaledResidual',norm(residual,inf), ...
    'massDefect_kgm',massDefect,'localMassDefect_kgm',cellMassDefect, ...
    'freeEnergyChange_Jm',deltaFreeEnergy,'chemicalWork_Jm',chemicalWork, ...
    'darcyDissipation_Jm',darcyDissipation,'phaseDissipation_Jm',phaseDissipation, ...
    'backwardEulerRemainder_Jm',convexRemainder,'energyIdentityDefect_Jm',energyIdentityDefect, ...
    'entropyProduced_JmK',entropyProduced,'thermostatHeat_Jm',thermostatHeat, ...
    'temperatureJumpHeat_Jm',temperatureJumpHeat,'temperatureJumpEntropy_JmK',temperatureJumpEntropy, ...
    'state',s,'flow',flow);
end

function update = newtonDirection(tangent,residual,z,scaling)
% Bound-constrained linearized least-squares correction. Retain ALL residual
% rows when fixing an outward phase direction, including the phase KKT row.
% The true KKT residual still determines line search and final convergence.
n = size(z,1); matrix = tangent*spdiags(scaling,0,3*n,3*n); rhs = -residual;
pinned = false(n,1); fixedTarget = zeros(n,1);
for attempt = 1:n+1
    fixed = 3*find(pinned); free = setdiff((1:3*n)',fixed);
    update = zeros(3*n,1); update(fixed) = fixedTarget(pinned);
    reducedRhs = rhs;
    if ~isempty(fixed), reducedRhs = rhs-matrix(:,fixed)*update(fixed); end
    freeMatrix = matrix(:,free);
    columnNorm = sqrt(full(sum(abs(freeMatrix).^2,1)))';
    assert(all(columnNorm>0),'JointInterface:SingularTangent','Zero tangent column detected.');
    columnScale = 1./columnNorm;
    balanced = freeMatrix*spdiags(columnScale,0,numel(free),numel(free));
    update(free) = columnScale.*(balanced\reducedRhs);
    phaseUpdate = update(3:3:end);
    lower = z(:,3)<=1e-13 & phaseUpdate < 0;
    upper = z(:,3)>=1-1e-13 & phaseUpdate > 0;
    add = (lower|upper) & ~pinned;
    if ~any(add), return; end
    pinned = pinned|add;
    fixedTarget(add) = double(upper(add))-z(add,3);
end
end

function settings = defaults(settings)
values = struct('length_m',0.45,'reservoirPressure_Pa',[2e5;0], ...
    'sealedEnds',[false;false],'outletSeepage',true,'maxIterations',60, ...
    'maxBacktracks',32,'tolerance',1e-10,'massAuditTolerance_kgm',1e-10,'trace',false, ...
    'energyAuditRelativeTolerance',1e-7);
names = fieldnames(values);
for k = 1:numel(names)
    if ~isfield(settings,names{k}), settings.(names{k}) = values.(names{k}); end
end
validateattributes(settings.length_m,{'double'},{'scalar','finite','positive'});
validateattributes(settings.reservoirPressure_Pa,{'double'},{'numel',2,'finite','real'});
assert(numel(settings.sealedEnds)==2,'JointInterface:BoundaryShape','Two end flags are required.');
settings.reservoirPressure_Pa = settings.reservoirPressure_Pa(:);
settings.sealedEnds = logical(settings.sealedEnds(:));
end

function alpha = feasibleStep(x,update,p)
alpha = 1;
constraints = [x(:,3),1-x(:,3),p.aperture0_m+p.wetAreaFraction*x(:,1)];
directions = [update(:,3),-update(:,3),p.wetAreaFraction*update(:,1)];
roundoff = 1e-14;
mask = directions < -roundoff;
if any(mask,'all'), alpha = min(alpha,min(-constraints(mask)./directions(mask))); end
alpha = max(0,alpha);
end

function ice = roundoffBounds(ice)
% Only floating-point roundoff at the KKT bounds is corrected, never pressure.
ice(ice<0 & ice>=-1e-13) = 0;
ice(ice>1 & ice<=1+1e-13) = 1;
end

function ok = admissible(x,temperature,p)
ok = temperature>0 && all(x(:,3)>=0 & x(:,3)<=1 & ...
    p.aperture0_m+p.wetAreaFraction*x(:,1)>0 & ...
    x(:,2)>=p.minimumPressure_Pa & ...
    abs(x(:,2))<=p.maximumCompression*min(p.bulkWater_Pa,p.bulkIce_Pa));
end
