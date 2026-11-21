function out = solve_thermodynamic_joint(options)
%SOLVE_THERMODYNAMIC_JOINT Conservative, thermostated slit verification model.
% Independent research branch: not automatically substituted into v16 THM.
% Temperature is uniform and imposed, wall compliance is a local condensed
% spring, ice is immobile, and the remaining planar liquid layer carries flow.
% No geometry damage, gasket failure, or nucleation barrier is invented.
arguments
    options (1,1) struct = struct()
end
defaults = struct('parameters',struct(),'cells',16,'length_m',0.45, ...
    'duration_s',365*86400,'dt_s',86400,'initialTemperature_K',281.15, ...
    'temperatureFunction',@(t)277.15+8*cos(2*pi*t/(365*86400)), ...
    'reservoirPressure_Pa',[2e5;0],'sealedEnds',[false;false], ...
    'outletSeepage',true,'initialState',[],'verbose',false,'failureDiagnosticFile',"", ...
    'minimumDt_s',1,'maximumRejectedSteps',2000,'timeGrid_s',[]);
names = fieldnames(options);
for k = 1:numel(names)
    assert(isfield(defaults,names{k}),'JointInterface:UnknownOption','Unknown option %s.',names{k});
    defaults.(names{k}) = options.(names{k});
end
o = defaults; p = joint_thermodynamic_parameters(o.parameters);
validateattributes(o.cells,{'double'},{'scalar','integer','positive'});
validateattributes(o.duration_s,{'double'},{'scalar','positive','finite'});
validateattributes(o.dt_s,{'double'},{'scalar','positive','finite'});
if ~isempty(o.timeGrid_s)
    o.timeGrid_s = o.timeGrid_s(:);
    assert(o.timeGrid_s(1)==0 && o.timeGrid_s(end)==o.duration_s && ...
        all(diff(o.timeGrid_s)>0),'JointInterface:TimeGrid','Prescribed time grid must cover the full duration.');
end
settings = struct('length_m',o.length_m,'reservoirPressure_Pa',o.reservoirPressure_Pa, ...
    'sealedEnds',o.sealedEnds,'outletSeepage',o.outletSeepage);
n = o.cells;
if isempty(o.initialState)
    pressure = o.reservoirPressure_Pa(1)+(o.reservoirPressure_Pa(2)- ...
        o.reservoirPressure_Pa(1))*((1:n)'-0.5)/n;
    if all(o.sealedEnds), pressure(:) = 0; end
    gap = pressure*p.wetAreaFraction/p.normalStiffness_Pam+p.referenceGap_m;
    compressed = gap<0;
    gap(compressed) = (pressure(compressed)*p.wetAreaFraction+ ...
        p.normalStiffness_Pam*p.referenceGap_m)/(p.normalStiffness_Pam+p.contactStiffness_Pam);
    mass = p.rhoWater*(p.aperture0_m+p.wetAreaFraction*gap)./(1-pressure/p.bulkWater_Pa);
    x = [gap,mass,zeros(n,1),repmat(o.initialTemperature_K,n,1)];
else
    x = o.initialState;
    assert(isequal(size(x),[n,4]),'JointInterface:InitialState','Initial state must be cells-by-4.');
end
joint_helmholtz_state(x,p);
initialState = x;
capacity = ceil(o.duration_s/o.dt_s)+100;
if ~isempty(o.timeGrid_s), capacity = numel(o.timeGrid_s)+100; end
times = zeros(capacity+1,1); states = zeros(n,4,capacity+1); states(:,:,1) = x;
pressureHistory = zeros(n,capacity+1);
s = joint_helmholtz_state(x,p); pressureHistory(:,1) = s.pressure_Pa;
outletFlux = zeros(capacity,1); inletFlux = outletFlux;
audits = zeros(capacity,9); durations = zeros(capacity,1); temperatures = zeros(capacity,1);
step = 1; proposedDt = o.dt_s; rejected = 0; rejectedLog = cell(0,3);
ticId = tic;
while times(step)<o.duration_s
    dt = min(proposedDt,o.duration_s-times(step));
    if ~isempty(o.timeGrid_s), dt = o.timeGrid_s(step+1)-o.timeGrid_s(step); end
    time = times(step)+dt;
    temperature = o.temperatureFunction(time);
    try
        [next,info] = joint_thermodynamic_step(x,temperature,dt,settings,p);
    catch exception
        numericalFailure = ismember(exception.identifier, ...
            {'JointInterface:LineSearch','JointInterface:NonConvergence'});
        if numericalFailure && isempty(o.timeGrid_s) && dt/2>=o.minimumDt_s && rejected<o.maximumRejectedSteps
            rejected = rejected+1;
            rejectedLog(end+1,:) = {time,dt,exception.identifier}; %#ok<AGROW>
            proposedDt = dt/2;
            continue % Old mass, phase and gap remain exactly unchanged.
        end
        if strlength(string(o.failureDiagnosticFile))>0
            failure = struct('time_s',time,'dt_s',dt,'old',x,'temperature_K',temperature, ...
                'settings',settings,'parameters',p,'message',getReport(exception,'extended','hyperlinks','off'));
            save(o.failureDiagnosticFile,'failure','-v7');
        end
        rethrow(exception)
    end
    x = next;
    if step>capacity
        capacity = capacity+100;
        times(capacity+1,1) = 0; states(n,4,capacity+1) = 0;
        pressureHistory(n,capacity+1) = 0;
        outletFlux(capacity,1) = 0; inletFlux(capacity,1) = 0;
        audits(capacity,9) = 0; durations(capacity,1) = 0; temperatures(capacity,1) = 0;
    end
    times(step+1) = time; states(:,:,step+1) = x;
    pressureHistory(:,step+1) = info.state.pressure_Pa;
    outletFlux(step) = info.flow.referenceVolumeFlux_m2s(end);
    inletFlux(step) = info.flow.referenceVolumeFlux_m2s(1);
    durations(step) = dt; temperatures(step) = temperature;
    audits(step,:) = [info.scaledResidual,info.massDefect_kgm, ...
        max(abs(info.localMassDefect_kgm)),info.energyIdentityDefect_Jm, ...
        info.darcyDissipation_Jm,info.phaseDissipation_Jm, ...
        info.backwardEulerRemainder_Jm,info.entropyProduced_JmK,info.iterations];
    if o.verbose && (mod(step,100)==0 || time==o.duration_s)
        fprintf('Thermostated slit: %.3f days, Q %.6g L/m, max p %.3f MPa\n', ...
            time/86400,1000*sum(outletFlux(1:step).*durations(1:step)),max(info.state.pressure_Pa)/1e6);
    end
    proposedDt = min(o.dt_s,dt*1.5);
    step = step+1;
end
steps = step-1;
times = times(1:steps+1); states = states(:,:,1:steps+1);
pressureHistory = pressureHistory(:,1:steps+1);
outletFlux = outletFlux(1:steps); inletFlux = inletFlux(1:steps);
audits = audits(1:steps,:); durations = durations(1:steps); temperatures = temperatures(1:steps);
dx = o.length_m/n;
storageChange = dx*sum(x(:,2)-initialState(:,2));
integratedMassIn = p.rhoWater*sum(inletFlux.*durations);
integratedMassOut = p.rhoWater*sum(outletFlux.*durations);
out = struct('engine','ThermostatedJoint_v1','scope', ...
    'Ideal planar slit; prescribed uniform temperature; local elastic wall compliance; no irreversible damage', ...
    'parameters',p,'options',o,'initialState',initialState,'finalState',x, ...
    'time_s',times,'dt_s',durations,'temperature_K',temperatures, ...
    'stateHistory',states,'pressureHistory_Pa',pressureHistory, ...
    'inletReferenceFlux_m2s',inletFlux,'outletReferenceFlux_m2s',outletFlux, ...
    'outletReferenceVolume_Lm',1000*sum(outletFlux.*durations), ...
    'massLedger_kgm',[storageChange,integratedMassIn,integratedMassOut, ...
        storageChange-integratedMassIn+integratedMassOut], ...
    'audits',audits,'auditColumns',{{'scaledResidual','globalMassDefect_kgm', ...
    'localMassDefect_kgm','energyIdentityDefect_Jm','darcyDissipation_Jm', ...
    'phaseDissipation_Jm','backwardEulerRemainder_Jm','entropyProduced_JmK','NewtonIterations'}}, ...
    'rejectedSteps',rejected,'rejectedLog',{rejectedLog},'elapsed_s',toc(ticId));
end
