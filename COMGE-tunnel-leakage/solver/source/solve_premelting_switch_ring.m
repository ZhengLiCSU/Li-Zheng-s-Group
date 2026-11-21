function out = solve_premelting_switch_ring(model,p,f,options)
%SOLVE_PREMELTING_SWITCH_RING Complete nonuniform thermal cycles, with history.
% Verification forcing is a synthetic accelerated cycle, NOT calendar years
% or a calibrated tunnel environment. Mode-I cohesive damage is irreversible; no calibrated gasket fatigue is claimed.
arguments
    model (1,1) struct
    p (1,1) struct = joint_threephase_parameters()
    f (1,1) struct = premelting_switch_parameters()
    options (1,1) struct = struct()
end
o=struct('cycles',2,'periodDays',20,'dtDays',0.5,'initialPressure_Pa',1e5, ...
    'reservoirPressure_Pa',2e5,'temperatureAmplitude_K',16,'sealed',true, ...
    'minimumStepSeconds',0.001,'saveFile','','timeGrid_s',[],'verbose',true, ...
    'maximumDamageIncrement',0.025,'maximumPhaseIncrement',0.05,'maximumPressureIncrement_Pa',1e5, ...
    'boundarySchedule',[],'maximumActiveIcePressureIncrement_Pa',Inf,'stepController','growth');
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(o,names{k}),'ThermoInterface:Option','Unknown run option %s.',names{k});
    o.(names{k})=options.(names{k});
end
assert(ismember(o.stepController,{'growth','transition_predictor'}), ...
    'ThermoInterface:TimeControl','Unknown step controller.');
positive={'cycles','periodDays','dtDays','minimumStepSeconds','maximumPhaseIncrement','maximumPressureIncrement_Pa','maximumDamageIncrement'};
for k=1:numel(positive)
    value=o.(positive{k});
    assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value) && value>0, ...
        'ThermoInterface:TimeControl','Run duration and transition controls must be finite positive scalars.');
end
assert(o.cycles==floor(o.cycles) && isscalar(o.sealed) && isscalar(o.verbose) && ...
    isscalar(o.reservoirPressure_Pa) && isfinite(o.reservoirPressure_Pa) && o.reservoirPressure_Pa>=0 && ...
    isscalar(o.initialPressure_Pa) && isfinite(o.initialPressure_Pa) && ...
    isscalar(o.temperatureAmplitude_K) && isfinite(o.temperatureAmplitude_K) && o.temperatureAmplitude_K>=0, ...
    'PhaseSwitch:Input','Cycles must be integer and forcing controls must be finite scalars.');
assert(model.referenceTemperature_K-o.temperatureAmplitude_K>=240,'PhaseSwitch:Input','Synthetic forcing exceeds the supported temperature range.');
assert(isempty(o.timeGrid_s) || (isvector(o.timeGrid_s) && all(isfinite(o.timeGrid_s)) && all(diff(o.timeGrid_s)>0) && o.timeGrid_s(end)>=o.cycles*o.periodDays*86400), ...
    'ThermoInterface:TimeGrid','Accepted grid must strictly increase and cover the run.');
assert(o.cycles>0 && o.periodDays>0 && o.dtDays>0 && o.minimumStepSeconds>0 && ...
    o.maximumPhaseIncrement>0 && o.maximumPhaseIncrement<=1 && o.maximumPressureIncrement_Pa>0 && ...
    isscalar(o.maximumActiveIcePressureIncrement_Pa) && o.maximumActiveIcePressureIncrement_Pa>0, ...
    'ThermoInterface:TimeControl','Durations and transition limits must be positive; phase increment cannot exceed one.');
assert(isempty(o.saveFile) || (~isfile(o.saveFile) && ~isfile([o.saveFile '.failure.mat'])), ...
    'ThermoInterface:OutputExists','Preserve existing completed or failed evidence; use a fresh output path.');
endTime=o.cycles*o.periodDays*86400;
if ~isempty(o.boundarySchedule)
    schedule=o.boundarySchedule;
    required={'time_s';'outerTemperature_K';'innerTemperature_K';'reservoirPressure_Pa'};
    assert(isstruct(schedule) && isequal(sort(fieldnames(schedule)),sort(required)), ...
        'PhaseSwitch:Input','Boundary schedule requires time_s, two temperatures in K and reservoir gauge pressure in Pa.');
    tgrid=schedule.time_s(:);
    assert(numel(tgrid)>=2 && all(isfinite(tgrid)) && all(diff(tgrid)>0) && tgrid(1)==0 && tgrid(end)>=endTime, ...
        'PhaseSwitch:Input','Boundary schedule must increase from zero and cover the whole run.');
    for k=2:numel(required)
        value=schedule.(required{k});
        assert(isvector(value) && numel(value)==numel(tgrid) && isreal(value) && all(isfinite(value)), ...
            'PhaseSwitch:Input','Boundary series have incompatible shapes or nonfinite values.');
    end
    assert(all(schedule.outerTemperature_K>=240 & schedule.outerTemperature_K<=310) && ...
        all(schedule.innerTemperature_K>=240 & schedule.innerTemperature_K<=310) && all(schedule.reservoirPressure_Pa>=0), ...
        'PhaseSwitch:Input','Boundary temperatures or pressure exceed the supported domain.');
end
model=prepare_premelting_solver(model,p);
n=numel(model.weight); ns=numel(model.capacity); T0=model.referenceTemperature_K;
assert(o.initialPressure_Pa+p.atmosphericPressure_Pa>0,'ThreePhase:Domain','Positive absolute initial pressure required.');
assert(p.phaseViscosity_Jm2skg2==0,'ThreePhase:Equilibrium','This branch uses equilibrium phase transitions.');
z=[o.initialPressure_Pa/1e6*ones(n,1);ones(n,1);zeros(2*n,1);T0*ones(n,1);o.initialPressure_Pa/1e6*ones(n,1);T0*ones(ns,1)];
history=zeros(n,1); initial=premelting_switch_ring_state(z,model,p,f,history);
history=initial.cohesive.kappa_m; historyArchive=history; damage=initial.cohesive.damage;
acceptedState=initial;
surfaceEnergy=sum(model.weight.*initial.cohesive.surfaceEnergy_Jm2); cohesiveDefect=zeros(0,1);
t=0; proposal=o.dtDays*86400;
times=0; states=z; audit=zeros(0,11); flow=zeros(size(model.faces,1),0);
reject=zeros(0,3); rejectReason=cell(0,1); iterationCount=0; timer=tic;
try
    while t<endTime-1e-6
        if isempty(o.timeGrid_s)
            cycleEnd=(floor(t/(o.periodDays*86400)+1e-10)+1)*o.periodDays*86400;
            dt=min([proposal,endTime-t,cycleEnd-t]);
            if ~isempty(o.boundarySchedule)
                knots=o.boundarySchedule.time_s(o.boundarySchedule.time_s>t+1e-6);
                if ~isempty(knots), dt=min(dt,knots(1)-t); end
            end
        else
            candidates=o.timeGrid_s(o.timeGrid_s>t+1e-6);
            assert(~isempty(candidates),'ThermoInterface:TimeGrid','Supplied grid does not cover the run.');
            dt=candidates(1)-t;
        end
        [outerT,innerT,head]=premelting_boundary_at(t+dt,T0,o);
        bc=thermo_interface_boundary(model,outerT,innerT,head,o.sealed);
        try
            [candidate,a]=premelting_switch_ring_step(z,dt,model,bc,p,f,history);
        catch exception
            recoverable=ismember(exception.identifier,{'ThermoInterface:LineSearch','ThermoInterface:NonConvergence','Cohesive:EnergyResolution','ThermoInterface:EnergyAudit'});
            if recoverable && isempty(o.timeGrid_s) && dt/2>=o.minimumStepSeconds
                reject(end+1,:)=[t,dt,numel(times)]; %#ok<AGROW>
                rejectReason{end+1,1}=exception.identifier; %#ok<AGROW>
                proposal=dt/2; continue
            end
            rethrow(exception)
        end
        phaseChange=max(abs(candidate(n+(1:3*n))-z(n+(1:3*n))));
        pressureChange=1e6*max(abs(candidate([1:n,5*n+(1:n)])-z([1:n,5*n+(1:n)])));
        damageChange=max(abs(a.state.cohesive.damage-damage(:,end)));
        transition=premelting_transition_change(acceptedState,a.state);
        if damageChange>o.maximumDamageIncrement || phaseChange>o.maximumPhaseIncrement || pressureChange>o.maximumPressureIncrement_Pa || ...
                transition.activeIcePressure_Pa>o.maximumActiveIcePressureIncrement_Pa
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
        acceptedState=a.state;
        audit(end+1,:)=[max(abs(a.localMassDefect_kgm)), ...
            max(abs(a.localFluidEnergyDefect_Jm)),max(abs(a.localSolidEnergyDefect_Jm)), ...
            a.totalEnergyDefect_Jm,a.entropyProduced_JmK,min(a.faceEntropyRate_WmK), ...
            norm(a.mechanicalResidual_Nm,inf),a.energyAcceptanceBound_Jm,a.boundaryEnergy_Jm, ...
            a.pressureWork_Jm,a.iterations]; %#ok<AGROW>
        flow(:,end+1)=a.faceFlux(:,1); %#ok<AGROW>
        iterationCount=iterationCount+a.iterations;
        proposal=premelting_step_proposal(dt,phaseChange,pressureChange,damageChange, ...
            transition.activeIcePressure_Pa,o);
        if o.verbose && mod(numel(times)-1,40)==0
            fprintf('Premelting three-phase day %.4f / %.4f; accepted %d rejected %d; ice %.3f..%.3f\n', ...
                t/86400,endTime/86400,numel(times)-1,size(reject,1),min(z(2*n+(1:n))),max(z(2*n+(1:n))));
        end
    end
catch exception
    failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
        'time_s',t,'lastAcceptedState',z,'time',times,'states',states,'audit',audit,'reject',reject, ...
        'model',model,'parameters',p,'filmParameters',f,'options',o,'rejectionReasons',{rejectReason}, ...
        'massFlux_kgms',flow,'surfaceEnergy_Jm',surfaceEnergy,'cohesiveQuadratureDefect_Jm',cohesiveDefect);
    failure.lastState=premelting_switch_ring_state(z,model,p,f,history); failure.history=historyArchive; failure.damage=damage;
    if ~isempty(o.saveFile), save([o.saveFile '.failure.mat'],'failure','-v7'); end
    rethrow(exception)
end
out=struct('engine','PremeltingSwitchRing_v3_exact_solver_cache','model',model,'parameters',p,'filmParameters',f,'options',o, ...
    'time_s',times,'states',states,'audit',audit,'massFlux_kgms',flow,'rejected',reject, ...
    'history',historyArchive,'damage',damage,'surfaceEnergy_Jm',surfaceEnergy,'cohesiveQuadratureDefect_Jm',cohesiveDefect, ...
    'wallSeconds',toc(timer),'newtonIterations',iterationCount,'completed',true,'rejectionReasons',{rejectReason}, ...
    'auditColumns',{{'localMass_kgm','localFluidEnergy_Jm','localSolidEnergy_Jm', ...
        'totalEnergy_Jm','entropyProduced_JmK','minFaceEntropy_WmK','mechanicalResidual_Nm', ...
        'energyAcceptanceBound_Jm','boundaryEnergy_Jm','pressureWork_Jm','iterations'}});
if ~isempty(o.saveFile), save(o.saveFile,'out','-v7'); end
end
