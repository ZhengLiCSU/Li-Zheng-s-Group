function [s,a] = premelting_prescribed_temperature_step(old,T,dt,model,bc,p,f)
%PREMELTING_PRESCRIBED_TEMPERATURE_STEP Conserved inventory, Q4 equilibrium.
% T and solid temperatures are prescribed thermal CONTROL variables here.
% This is NOT a solution of the transient heat PDE. Required thermostat heat
% is returned as a reaction, not presented as an independent energy test.
% Previous accepted cohesive history is held fixed throughout Newton.
arguments
    old (1,1) struct
    T (:,1) double {mustBeFinite}
    dt (1,1) double {mustBeFinite,mustBePositive}
    model (1,1) struct
    bc (1,1) struct
    p (1,1) struct
    f (1,1) struct
end
pl=old.slit.equilibrium.liquidPressure_Pa; n=numel(pl); w=model.weight;
history=old.cohesive.kappa_m; Ts=old.solidTemperature_K;
for iteration=1:35
    [s,d]=premelting_ring_state(pl,T,Ts,model,p,f,history);
    [flux,JF]=premelting_ring_flux(s,model,bc,p);
    r=w.*(s.slit.mass_kgm2-old.slit.mass_kgm2)+dt*model.divergence*flux(:,1);
    % Full state representability includes absolute temperature, not just
    % stored mass: dM/dT can be large in a mechanically compliant slit.
    roundoff=8*w.*(abs(d.mass)*eps(abs([pl;T;Ts])));
    bound=max(roundoff,64*w.*eps(max(abs(s.slit.mass_kgm2),abs(old.slit.mass_kgm2))));
    if all(abs(r)<=bound), break; end
    J=w.*d.mass(:,1:n)+dt*model.divergence*JF;
    update=-(J\r); alpha=1; accepted=false;
    for k=1:35
        try
            trial=premelting_ring_state(pl+alpha*update,T,Ts,model,p,f,history);
            fl=premelting_ring_flux(trial,model,bc,p);
            rt=w.*(trial.slit.mass_kgm2-old.slit.mass_kgm2)+dt*model.divergence*fl(:,1);
            valid=true;
        catch exception
            if startsWith(exception.identifier,'Premelting:')||startsWith(exception.identifier,'ThreePhase:')
                valid=false;
            else
                rethrow(exception)
            end
        end
        if valid && (norm(rt,inf)<(1-1e-4*alpha)*norm(r,inf)||all(abs(rt)<=bound))
            pl=pl+alpha*update; accepted=true; break
        end
        alpha=alpha/2;
    end
    assert(accepted,'Premelting:InventoryDomain', ...
        'No admissible conserved-inventory update. Cavitation/topology change may need the full three-phase model.');
end
assert(all(abs(r)<=bound),'Premelting:NonConvergence','Conserved inventory solve failed.');
assert(max(abs(r))<1e-12,'Premelting:MassAudit','Independent physical mass gate failed.');
[flux,~,ports]=premelting_ring_flux(s,model,bc,p);
deltaGap=s.gap_m-old.gap_m; piAvg=(s.slit.equilibrium.icePressure_Pa+old.slit.equilibrium.icePressure_Pa)/2;
work=w*p.wetAreaFraction.*piAvg.*deltaGap;
heat=w.*(s.slit.internalEnergy_Jm2-old.slit.internalEnergy_Jm2)+work+dt*model.divergence*flux(:,2);
du=s.u-old.u;
solidEnergy=0.5*(s.u+old.u)'*model.stiffness*du+model.referenceTemperature_K*sum(model.thermalLoad'*du)+ ...
    sum(w.*(s.cohesive.elasticEnergy_Jm2-old.cohesive.elasticEnergy_Jm2+ ...
    s.cohesive.surfaceEnergy_Jm2-old.cohesive.surfaceEnergy_Jm2));
% Holding the host temperature fixed also needs reciprocal thermoelastic heat.
hostHeat=sum(Ts.*(model.thermalLoad'*du));
mechanicalWorkDefect=solidEnergy-hostHeat-sum(work);
a=struct('localMassDefect_kgm',r,'localMassAcceptanceBound_kgm',bound, ...
    'faceFlux',flux,'faceEntropyRate_WmK',ports.faceEntropyRate_WmK, ...
    'thermostatHeat_Jm',heat,'hostThermostatHeat_Jm',hostHeat, ...
    'mechanicalWorkQuadratureDefect_Jm',mechanicalWorkDefect,'iterations',iteration, ...
    'scope','Prescribed-temperature inventory/FE solve; thermostat heat is a reaction, not a solved heat balance.');
end
