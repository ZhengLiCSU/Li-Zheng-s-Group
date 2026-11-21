function q = summarize_premelting_sample(out)
%SUMMARIZE_PREMELTING_SAMPLE Integrate accepted fluxes, preserving eight joints.
assert(out.completed,'PhaseSwitch:Incomplete','Only a complete trajectory has a leakage QoI.');
m=out.model; p=out.parameters; dt=diff(out.time_s);
mult=1; if isfield(m,'symmetry'), mult=m.symmetry.multiplicity; end
outlet=m.faces(:,2)==0; inlet=m.faces(:,1)==0;
jointQ=sum(out.massFlux_kgms(outlet,:).*dt',2)/p.rhoWater*1000;
if mult>1, jointQ=repmat(jointQ,mult,1); end
assert(numel(jointQ)==8,'PhaseSwitch:Shape','Eight longitudinal joints are required.');
cycleEnd=out.options.cycles*out.options.periodDays*86400;
start=cycleEnd-out.options.periodDays*86400;
overlap=max(0,min(out.time_s(2:end),cycleEnd)-max(out.time_s(1:end-1),start));
last=sum(out.massFlux_kgms(outlet,:).*overlap',2)/p.rhoWater*1000;
if mult>1, last=repmat(last,mult,1); end
ice=zeros(numel(out.time_s),1); vapor=ice; peak=-Inf; minT=Inf; maxStrain=0;
for k=1:numel(out.time_s)
    s=premelting_switch_ring_state(out.states(:,k),m,p,out.filmParameters,out.history(:,k));
    ice(k)=mult*sum(m.weight.*s.iceVolume_m); vapor(k)=mult*sum(m.weight.*s.vaporVolume_m);
    active=s.areaFraction(:,2)>1e-10;
    if any(active), peak=max(peak,max(s.eos.pressure(active,2),[],'all')); end
    minT=min(minT,min(s.x(:,4)));
    strain=reshape(m.strainOperator*s.u,3,[])';
    maxStrain=max(maxStrain,max(abs([strain(:,1),strain(:,2),strain(:,3)]),[],'all'));
end
if ~isfinite(peak), peak=NaN; end % Absence of ice is not zero ice pressure.
q=struct('caseCompleted',true,'independentAuditPending',true, ...
    'leakage_Lpm',sum(jointQ),'jointLeakage_Lpm',jointQ, ...
    'lastCycleLeakage_Lpm',sum(last),'lastCycleJointLeakage_Lpm',last, ...
    'inflow_Lpm',mult*sum(out.massFlux_kgms(inlet,:).*dt','all')/p.rhoWater*1000, ...
    'peakActiveIcePressure_Pa',peak,'minimumJointTemperature_K',minT, ...
    'timeIntegratedIceVolume_m2day',sum(ice(2:end).*dt)/86400, ...
    'timeIntegratedVaporVolume_m2day',sum(vapor(2:end).*dt)/86400, ...
    'maximumStrainComponent',maxStrain,'acceptedSteps',numel(dt), ...
    'rejectedSteps',size(out.rejected,1),'wallSeconds',out.wallSeconds, ...
    'symmetryMultiplicity',mult,'jointBodyLeakage_Lpm',0, ...
    'scope','Component cycle, L per metre of tunnel axis; not a calibrated calendar-year prediction.');
end
