function [next,audit] = premelting_switch_ring_step(old,dt,model,bc,p,f,history,options)
%PREMELTING_SWITCH_RING_STEP Semismooth simplex Newton with physical audits.
% All accepted states obey pressure/phase/temperature domains, local mass,
% local heat, global energy and a POSTERIORI global entropy audit.
arguments
    old (:,1) double {mustBeFinite}
    dt (1,1) double {mustBePositive,mustBeFinite}
    model (1,1) struct
    bc (1,1) struct
    p (1,1) struct
    f (1,1) struct
    history (:,1) double {mustBeFinite,mustBeNonnegative}
    options (1,1) struct = struct()
end
v=struct('tolerance',1e-11,'maximumIterations',60,'trace',false,'eliminateTrivialPhaseRows',true);
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(v,names{k}),'ThermoInterface:Option','Unknown solver option.');
    v.(names{k})=options.(names{k});
end
model=prepare_premelting_solver(model,p);
n=numel(model.weight); ns=numel(model.capacity); nz=numel(old);
scaling=ones(nz,1);
next=old; iterations=0;
oldState=premelting_switch_ring_state(old,model,p,f,history);
for iteration=1:v.maximumIterations
    [r,J,current]=premelting_switch_ring_residual(next,old,dt,model,bc,p,f,history,[],oldState);
    % A representability bound, derived from the CURRENT exact Jacobian and
    % one ulp of each absolute variable, prevents futile substeps when the
    % requested residual lies below temperature floating-point resolution.
    % Conservation audits below remain independent and unchanged.
    residualBound=max(v.tolerance,4*abs(J)*eps(abs(next)));
    % A componentwise scaled residual can pass near steady state while the
    % sum of small heat residuals still violates the physical global ledger.
    % Continue Newton in that case; NEVER relax the conservation threshold.
    sid=6*n+(1:ns);
    floorE=16*sum(model.capacity.*eps(max(abs(next(sid)),abs(old(sid)))));
    scaleE=max([1,abs(current.boundaryEnergy_Jm),abs(current.solidEnergyChange_Jm),abs(current.fluidEnergyChange_Jm)]);
    closed=abs(current.totalEnergyDefect_Jm)<max(1e-8*scaleE,floorE) && ...
        max(abs(current.localMassDefect_kgm))<1e-10 && ...
        max(abs(current.localFluidEnergyDefect_Jm))<1e-5 && ...
        max(abs(current.localSolidEnergyDefect_Jm))<1e-3;
    if all(abs(r)<=residualBound) && closed, break; end
    if all(abs(r)<=residualBound) && abs(current.cohesiveQuadratureDefect_Jm)>max(1e-8*scaleE,floorE)/2
        error('Cohesive:EnergyResolution','Resolve a cohesive loading/unloading kink by smaller time steps.');
    end
    if v.trace, fprintf('Coupled ring Newton %d R=%.4g\n',iteration,norm(r,inf)); end
    loading=[]; consistent=false;
    for branchIteration=1:3*n+2
        if ~isempty(loading)
            [r,J]=premelting_switch_ring_residual(next,old,dt,model,bc,p,f,history,loading,oldState);
        end
        matrix=J.*scaling'; rhs=-r;
        update=premelting_newton_direction(matrix,rhs,next(n+(1:3*n)),n,v.eliminateTrivialPhaseRows);
        % An elastic interface has no history-tangent kink to classify.
        if ~model.cohesive.enabled, consistent=true; break; end
        [state,derivative]=premelting_switch_ring_state(next,model,p,f,history,loading);
        ties=cohesive_history_ties(state.x(:,1),history,model.cohesive);
        gapDirection=derivative.gap*(scaling.*update);
        proposedLoading=gapDirection>=0;
        currentLoading=state.cohesive.tangent_Pam<0;
        % A fully frozen, hydraulically isolated cell can have exactly zero
        % gap increment. Do not toggle its tangent on a roundoff-only sign.
        directionFloor=4*abs(derivative.gap)*eps(abs(next));
        unresolved=abs(gapDirection)<=directionFloor;
        proposedLoading(unresolved)=currentLoading(unresolved);
        if all(proposedLoading(ties)==currentLoading(ties)), consistent=true; break; end
        loading=proposedLoading;
    end
    assert(consistent,'ThermoInterface:NonConvergence','Loading/unloading directional active set did not stabilize.');
    update=scaling.*update; alpha=1;
    fraction=next(n+(1:3*n)); delta=update(n+(1:3*n));
    if any(delta<0), alpha=min(alpha,min(-fraction(delta<0)./delta(delta<0))); end
    if any(delta>0), alpha=min(alpha,min((1-fraction(delta>0))./delta(delta>0))); end
    accepted=false;
    for backtrack=1:35
        trial=next+alpha*update; phase=trial(n+(1:3*n));
        phase(phase<0 & phase>=-1e-13)=0; phase(phase>1 & phase<=1+1e-13)=1;
        trial(n+(1:3*n))=phase;
        try
            rt=premelting_switch_ring_residual(trial,old,dt,model,bc,p,f,history,[],oldState);
            valid=true;
        catch exception
            if startsWith(exception.identifier,'PhaseSwitch:') || startsWith(exception.identifier,'ThreePhase:') || startsWith(exception.identifier,'JointThermo:') || startsWith(exception.identifier,'WaterViscosity:') || strcmp(exception.identifier,'JointNonisothermal:Temperature') || strcmp(exception.identifier,'Cohesive:Mechanics')
                valid=false;
            else
                rethrow(exception)
            end
        end
        if v.trace
            if valid, fprintf('  line %d alpha %.3e R %.4e\n',backtrack,alpha,norm(rt,2));
            else, fprintf('  line %d alpha %.3e invalid: %s\n',backtrack,alpha,exception.message); end
        end
        if valid && (norm(rt,2)<=(1-1e-4*alpha)*norm(r,2) || all(abs(rt)<=residualBound))
            next=trial; accepted=true; break
        end
        alpha=alpha/2;
    end
    assert(accepted,'ThermoInterface:LineSearch','No feasible Newton step; R=%.4g.',norm(r,inf));
    iterations=iteration;
end
[r,J,audit]=premelting_switch_ring_residual(next,old,dt,model,bc,p,f,history,[],oldState);
residualBound=max(v.tolerance,4*abs(J)*eps(abs(next)));
assert(all(abs(r)<=residualBound),'ThermoInterface:NonConvergence','Residual %.4g exceeds its verified precision envelope.',norm(r,inf));
audit.maximumResidualRoundoffBound=max(residualBound);
if isfield(model,'strainOperator')
    strain=reshape(model.strainOperator*audit.state.u,3,[])';
    principalMean=(strain(:,1)+strain(:,2))/2;
    radius=hypot((strain(:,1)-strain(:,2))/2,strain(:,3)/2);
    audit.maximumPrincipalStrain=max(abs([principalMean+radius;principalMean-radius]));
    assert(audit.maximumPrincipalStrain<=0.05,'ThermoInterface:StrainDomain', ...
        'Small-strain host domain exceeded; finite-deformation mechanics is required.');
end
assert(max(abs(audit.localMassDefect_kgm))<1e-10,'ThermoInterface:MassAudit','Local mass audit failed.');
assert(max(abs(audit.localFluidEnergyDefect_Jm))<1e-5 && ...
    max(abs(audit.localSolidEnergyDefect_Jm))<1e-3,'ThermoInterface:LocalEnergyAudit','Local energy audit failed.');
energyScale=max([1,abs(audit.boundaryEnergy_Jm),abs(audit.solidEnergyChange_Jm),abs(audit.fluidEnergyChange_Jm)]);
% Subtracting two absolute temperatures in a massive solid has a known
% floating-point energy floor. Report it separately; NEVER zero the defect.
solidIDs=6*n+(1:ns);
roundoffBound=16*sum(model.capacity.*eps(max(abs(next(solidIDs)),abs(old(solidIDs)))));
audit.energyRoundoffBound_Jm=roundoffBound;
audit.energyAcceptanceBound_Jm=max(1e-8*energyScale,roundoffBound);
assert(abs(audit.totalEnergyDefect_Jm)<audit.energyAcceptanceBound_Jm, ...
    'ThermoInterface:EnergyAudit','Total energy defect %.5g exceeds %.5g J/m.', ...
    audit.totalEnergyDefect_Jm,audit.energyAcceptanceBound_Jm);
assert(audit.entropyProduced_JmK>=-1e-7 && min(audit.faceEntropyRate_WmK)>=-1e-12, ...
    'ThermoInterface:EntropyAudit','Negative discrete entropy production; refine or diagnose, do not clip.');
assert(all(audit.state.cohesive.kappa_m>=history), ...
    'Cohesive:Healing','Irreversible history decreased.');
assert(audit.fractureSurfaceEnergyChange_Jm>=-1e-12,'Cohesive:Healing','Fracture surface energy decreased.');
assert(norm(audit.mechanicalResidual_Nm,inf)<1e-4,'Cohesive:MechanicalAudit','Nonlinear mechanical equilibrium failed.');
assert(max(abs(sum(audit.state.areaFraction,2)-1))<1e-10 && ...
    all(audit.state.areaFraction>=-1e-13 & audit.state.areaFraction<=1+1e-13,'all'), ...
    'PhaseSwitch:Simplex','Accepted phase areas do not form a simplex.');
activeIce=audit.state.areaFraction(:,2)>1e-10;
assert(all(audit.state.filmBranch(activeIce)==0) && ...
    all(2*audit.state.film_m(activeIce)./audit.state.aperture_m(activeIce)<f.maximumFilmFraction), ...
    'PhaseSwitch:FilmDomain','Active ice left the separated continuum-film domain; no clipping permitted.');
assert(max(abs(audit.wallResidual_Pa))<1e-3,'PhaseSwitch:WallAudit','Wall force conjugacy failed.');
audit.iterations=iterations;
end
