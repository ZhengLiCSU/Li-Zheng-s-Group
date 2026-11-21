function [next,audit] = thermo_interface_step(old,dt,model,bc,p,options)
%THERMO_INTERFACE_STEP Monolithic Newton solve with exact block derivatives.
% All accepted states obey pressure/phase/temperature domains, local mass,
% local heat, global energy and a POSTERIORI global entropy audit.
arguments
    old (:,1) double {mustBeFinite}
    dt (1,1) double {mustBePositive,mustBeFinite}
    model (1,1) struct
    bc (1,1) struct
    p (1,1) struct = joint_thermodynamic_parameters()
    options (1,1) struct = struct()
end
v=struct('tolerance',1e-11,'maximumIterations',60,'trace',false);
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(v,names{k}),'ThermoInterface:Option','Unknown solver option.');
    v.(names{k})=options.(names{k});
end
n=numel(model.weight); ns=numel(model.capacity); nz=numel(old);
scaling=[1e6*ones(n,1);ones(n,1);ones(n+ns,1)];
next=old; iterations=0;
for iteration=1:v.maximumIterations
    [r,J]=thermo_interface_residual(next,old,dt,model,bc,p);
    % A representability bound, derived from the CURRENT exact Jacobian and
    % one ulp of each absolute variable, prevents futile substeps when the
    % requested residual lies below temperature floating-point resolution.
    % Conservation audits below remain independent and unchanged.
    residualBound=max(v.tolerance,4*abs(J)*eps(abs(next)));
    if all(abs(r)<=residualBound), break; end
    if v.trace, fprintf('Coupled ring Newton %d R=%.4g\n',iteration,norm(r,inf)); end
    matrix=J.*scaling'; rhs=-r; pinned=false(n,1); target=zeros(n,1);
    for attempt=1:n+1
        fixed=n+find(pinned); free=setdiff((1:nz)',fixed);
        update=zeros(nz,1); update(fixed)=target(pinned);
        reduced=rhs;
        if ~isempty(fixed), reduced=rhs-matrix(:,fixed)*update(fixed); end
        columnNorm=sqrt(sum(matrix(:,free).^2,1))';
        update(free)=(matrix(:,free)./columnNorm'\reduced)./columnNorm;
        phaseUpdate=update(n+(1:n)); fraction=next(n+(1:n));
        lower=fraction<=1e-13 & phaseUpdate<0; upper=fraction>=1-1e-13 & phaseUpdate>0;
        add=(lower|upper)&~pinned;
        if ~any(add), break; end
        pinned=pinned|add; target(add)=double(upper(add))-fraction(add);
    end
    update=scaling.*update; alpha=1;
    fraction=next(n+(1:n)); delta=update(n+(1:n));
    if any(delta<0), alpha=min(alpha,min(-fraction(delta<0)./delta(delta<0))); end
    if any(delta>0), alpha=min(alpha,min((1-fraction(delta>0))./delta(delta>0))); end
    accepted=false;
    for backtrack=1:35
        trial=next+alpha*update; phase=trial(n+(1:n));
        phase(phase<0 & phase>=-1e-13)=0; phase(phase>1 & phase<=1+1e-13)=1;
        trial(n+(1:n))=phase;
        try
            rt=thermo_interface_residual(trial,old,dt,model,bc,p);
            valid=true;
        catch exception
            if startsWith(exception.identifier,'JointThermo:') || strcmp(exception.identifier,'JointNonisothermal:Temperature')
                valid=false;
            else
                rethrow(exception)
            end
        end
        if valid && (norm(rt,2)<=(1-1e-4*alpha)*norm(r,2) || all(abs(rt)<=residualBound))
            next=trial; accepted=true; break
        end
        alpha=alpha/2;
    end
    assert(accepted,'ThermoInterface:LineSearch','No feasible Newton step; R=%.4g.',norm(r,inf));
    iterations=iteration;
end
[r,J,audit]=thermo_interface_residual(next,old,dt,model,bc,p);
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
solidIDs=3*n+(1:ns);
roundoffBound=16*sum(model.capacity.*eps(max(abs(next(solidIDs)),abs(old(solidIDs)))));
audit.energyRoundoffBound_Jm=roundoffBound;
audit.energyAcceptanceBound_Jm=max(1e-8*energyScale,roundoffBound);
assert(abs(audit.totalEnergyDefect_Jm)<audit.energyAcceptanceBound_Jm, ...
    'ThermoInterface:EnergyAudit','Total energy defect %.5g exceeds %.5g J/m.', ...
    audit.totalEnergyDefect_Jm,audit.energyAcceptanceBound_Jm);
assert(audit.entropyProduced_JmK>=-1e-7 && min(audit.faceEntropyRate_WmK)>=-1e-12, ...
    'ThermoInterface:EntropyAudit','Negative discrete entropy production; refine or diagnose, do not clip.');
audit.iterations=iterations;
end
