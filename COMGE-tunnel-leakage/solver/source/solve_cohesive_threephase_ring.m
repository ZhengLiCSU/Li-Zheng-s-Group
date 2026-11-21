function out = solve_cohesive_threephase_ring(model,p,options)
%SOLVE_NONISOTHERMAL_RING Complete nonuniform thermal cycles, with history.
% Verification forcing is a synthetic accelerated cycle, NOT calendar years
% or a calibrated tunnel environment. Mode-I cohesive damage is irreversible; no calibrated gasket fatigue is claimed.
arguments
    model (1,1) struct
    p (1,1) struct = joint_threephase_parameters()
    options (1,1) struct = struct()
end
o=struct('cycles',2,'periodDays',20,'dtDays',0.125,'initialPressure_Pa',1e5, ...
    'reservoirPressure_Pa',2e5,'temperatureAmplitude_K',10,'sealed',false, ...
    'minimumStepSeconds',0.001,'saveFile','','timeGrid_s',[],'verbose',true, ...
    'maximumDamageIncrement',0.025,'maximumLogPressureIncrement',0.25,'maximumPhaseIncrement',0.025,'maximumPressureIncrement_Pa',1e4);
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(o,names{k}),'ThermoInterface:Option','Unknown run option %s.',names{k});
    o.(names{k})=options.(names{k});
end
assert(o.cycles>0 && o.periodDays>0 && o.dtDays>0 && o.minimumStepSeconds>0 && ...
    o.maximumPhaseIncrement>0 && o.maximumPhaseIncrement<=1 && o.maximumPressureIncrement_Pa>0, ...
    'ThermoInterface:TimeControl','Durations and transition limits must be positive; phase increment cannot exceed one.');
assert(isempty(o.saveFile) || (~isfile(o.saveFile) && ~isfile([o.saveFile '.failure.mat'])), ...
    'ThermoInterface:OutputExists','Preserve existing completed or failed evidence; use a fresh output path.');
n=numel(model.weight); ns=numel(model.capacity); T0=model.referenceTemperature_K;
assert(o.initialPressure_Pa+p.atmosphericPressure_Pa>0,'ThreePhase:Domain','Positive absolute initial pressure required.');
assert(p.phaseViscosity_Jm2skg2==0,'ThreePhase:Equilibrium','This branch uses equilibrium phase transitions.');
z=[log1p(o.initialPressure_Pa/p.atmosphericPressure_Pa)*ones(n,1);zeros(2*n,1);T0*ones(n+ns,1)];
history=zeros(n,1); initial=cohesive_threephase_state(z,model,p,history);
history=initial.cohesive.kappa_m; historyArchive=history; damage=initial.cohesive.damage;
surfaceEnergy=sum(model.weight.*initial.cohesive.surfaceEnergy_Jm2); cohesiveDefect=zeros(0,1);
t=0; endTime=o.cycles*o.periodDays*86400; proposal=o.dtDays*86400;
times=0; states=z; audit=zeros(0,11); flow=zeros(size(model.faces,1),0);
reject=zeros(0,3); rejectReason=cell(0,1); iterationCount=0; timer=tic;
try
    while t<endTime-1e-6
        if isempty(o.timeGrid_s)
            cycleEnd=(floor(t/(o.periodDays*86400)+1e-10)+1)*o.periodDays*86400;
            dt=min([proposal,endTime-t,cycleEnd-t]);
        else
            candidates=o.timeGrid_s(o.timeGrid_s>t+1e-6);
            assert(~isempty(candidates),'ThermoInterface:TimeGrid','Supplied grid does not cover the run.');
            dt=candidates(1)-t;
        end
        phase=(t+dt)/(o.periodDays*86400);
        innerT=T0-o.temperatureAmplitude_K*sin(pi*phase)^2;
        bc=thermo_interface_boundary(model,T0,innerT,o.reservoirPressure_Pa,o.sealed);
        try
            [candidate,a]=cohesive_threephase_step(z,dt,model,bc,p,history);
        catch exception
            recoverable=ismember(exception.identifier,{'ThermoInterface:LineSearch','ThermoInterface:NonConvergence','Cohesive:EnergyResolution','ThermoInterface:EnergyAudit'});
            if recoverable && isempty(o.timeGrid_s) && dt/2>=o.minimumStepSeconds
                reject(end+1,:)=[t,dt,numel(times)]; %#ok<AGROW>
                rejectReason{end+1,1}=exception.identifier; %#ok<AGROW>
                proposal=dt/2; continue
            end
            rethrow(exception)
        end
        phaseChange=max(abs(candidate(n+(1:2*n))-z(n+(1:2*n))));
        pressureChange=p.atmosphericPressure_Pa*max(abs(exp(candidate(1:n))-exp(z(1:n))));
        logPressureChange=max(abs(candidate(1:n)-z(1:n)));
        damageChange=max(abs(a.state.cohesive.damage-damage(:,end)));
        if damageChange>o.maximumDamageIncrement || phaseChange>o.maximumPhaseIncrement || pressureChange>o.maximumPressureIncrement_Pa || logPressureChange>o.maximumLogPressureIncrement
            assert(isempty(o.timeGrid_s) && dt/2>=o.minimumStepSeconds, ...
                'ThermoInterface:TransitionResolution','Phase/pressure increment needs a finer accepted time grid.');
            reject(end+1,:)=[t,dt,numel(times)]; %#ok<AGROW>
            rejectReason{end+1,1}='ThermoInterface:TransitionResolution'; %#ok<AGROW>
            proposal=dt/2; continue
        end
        history=a.state.cohesive.kappa_m; historyArchive(:,end+1)=history; %#ok<AGROW>
        damage(:,end+1)=a.state.cohesive.damage; surfaceEnergy(end+1,1)=a.surfaceEnergy_Jm; %#ok<AGROW>
        cohesiveDefect(end+1,1)=a.cohesiveQuadratureDefect_Jm; %#ok<AGROW>
        z=candidate; t=t+dt; times(end+1,1)=t; states(:,end+1)=z; %#ok<AGROW>
        audit(end+1,:)=[max(abs(a.localMassDefect_kgm)), ...
            max(abs(a.localFluidEnergyDefect_Jm)),max(abs(a.localSolidEnergyDefect_Jm)), ...
            a.totalEnergyDefect_Jm,a.entropyProduced_JmK,min(a.faceEntropyRate_WmK), ...
            norm(a.mechanicalResidual_Nm,inf),a.energyAcceptanceBound_Jm,a.boundaryEnergy_Jm, ...
            a.pressureWork_Jm,a.iterations]; %#ok<AGROW>
        flow(:,end+1)=a.faceFlux(:,1); %#ok<AGROW>
        iterationCount=iterationCount+a.iterations;
        proposal=min(o.dtDays*86400,1.5*dt);
        if o.verbose && mod(numel(times)-1,40)==0
            fprintf('Nonisothermal ring day %.4f / %.4f; accepted %d rejected %d; ice %.3f..%.3f\n', ...
                t/86400,endTime/86400,numel(times)-1,size(reject,1),min(z(n+(1:n))),max(z(n+(1:n))));
        end
    end
catch exception
    failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
        'time_s',t,'lastAcceptedState',z,'time',times,'states',states,'audit',audit,'reject',reject, ...
        'model',model,'parameters',p,'options',o,'rejectionReasons',{rejectReason});
    failure.lastState=cohesive_threephase_state(z,model,p,history); failure.history=historyArchive; failure.damage=damage;
    if ~isempty(o.saveFile), save([o.saveFile '.failure.mat'],'failure','-v7'); end
    rethrow(exception)
end
out=struct('engine','CohesiveThreePhaseRing_v1_surface_energy','model',model,'parameters',p,'options',o, ...
    'time_s',times,'states',states,'audit',audit,'massFlux_kgms',flow,'rejected',reject, ...
    'history',historyArchive,'damage',damage,'surfaceEnergy_Jm',surfaceEnergy,'cohesiveQuadratureDefect_Jm',cohesiveDefect, ...
    'wallSeconds',toc(timer),'newtonIterations',iterationCount,'completed',true,'rejectionReasons',{rejectReason}, ...
    'auditColumns',{{'localMass_kgm','localFluidEnergy_Jm','localSolidEnergy_Jm', ...
        'totalEnergy_Jm','entropyProduced_JmK','minFaceEntropy_WmK','mechanicalResidual_Nm', ...
        'energyAcceptanceBound_Jm','boundaryEnergy_Jm','pressureWork_Jm','iterations'}});
if ~isempty(o.saveFile), save(o.saveFile,'out','-v7'); end
end
