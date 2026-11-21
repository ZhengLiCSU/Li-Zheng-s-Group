function report = run_threephase_ring_verification(folder)
%RUN_THREEPHASE_RING_VERIFICATION Reproduce the limited three-phase campaign.
% New directory REQUIRED. Two synthetic 20-day cycles, not two years.
% Independent reconstruction/figures: audit_threephase_ring.py and
% plot_threephase_ring.py. No legacy 20-year source or results are changed.
arguments
    folder (1,:) char
end
assert(~isfolder(folder) && ~isfile(folder),'ThreePhase:OutputExists', ...
    'Use a new folder; do not overwrite completed or failed evidence.');
mkdir(folder); snapshot=fullfile(folder,'source_at_start'); mkdir(snapshot);
source=fileparts(mfilename('fullpath')); files=dir(fullfile(source,'*.m'));
for k=1:numel(files)
    copyfile(fullfile(source,files(k).name),fullfile(snapshot,files(k).name));
end
report=struct('completed',false,'scope','Pure water, equilibrium planar phases, no air or damage; synthetic cycles.');
tests=runtests(fullfile(source,'TestThreePhaseInterface.m'));
save(fullfile(folder,'threephase_tests.mat'),'tests');
assert(all([tests.Passed]),'ThreePhase:Tests','Material and coupling tests must pass.');
p=joint_threephase_parameters();
m=build_thermo_interface_ring(make_thm_verification_mesh(16),3,p);
coarse=solve_threephase_ring(m,p,struct('temperatureAmplitude_K',20, ...
    'saveFile',fullfile(folder,'strong_thermal_reference.mat')));
fineOptions=struct('temperatureAmplitude_K',20,'dtDays',0.0625, ...
    'maximumPhaseIncrement',0.0125,'maximumPressureIncrement_Pa',5000, ...
    'maximumLogPressureIncrement',0.125,'saveFile',fullfile(folder,'strong_thermal_fine.mat'));
fine=solve_threephase_ring(m,p,fineOptions);
fineOptions.temperatureAmplitude_K=0;
fineOptions.saveFile=fullfile(folder,'warm_control_energy_gate.mat');
warm=solve_threephase_ring(m,p,fineOptions);
report.reference=metrics(coarse); report.fine=metrics(fine); report.warm=metrics(warm);
a=report.reference; b=report.fine;
report.relativeErrors=[abs(a.cycleVolume-b.cycleVolume)/b.cycleVolume, ...
    abs(a.peakPressure-b.peakPressure)/b.peakPressure, ...
    abs(a.peakGap-b.peakGap)/b.peakGap, ...
    abs(a.peakVapor-b.peakVapor)/b.peakVapor, ...
    max(abs(a.bins-b.bins))/b.cycleVolume];
report.criteriaPassed=all(report.relativeErrors<=0.01);
report.completed=true;
save(fullfile(folder,'campaign.mat'),'report');
assert(report.criteriaPassed,'ThreePhase:TimeResolution','Preserve evidence: target-cycle time refinement exceeds 1%%.');
end

function a=metrics(out)
n=numel(out.model.weight); p=out.parameters; t=out.time_s(:)';
period=out.options.periodDays*86400;
outlet=out.model.faces(:,2)==0;
rate=sum(out.massFlux_kgms(outlet,:),1)*1000/p.rhoWater;
bins=linspace(period,2*period,13); volume=zeros(1,12);
for k=1:12
    duration=max(0,min(t(2:end),bins(k+1))-max(t(1:end-1),bins(k)));
    volume(k)=sum(rate.*duration);
end
selected=t>=period-1e-6; z=out.states(:,selected);
pg=p.atmosphericPressure_Pa*expm1(z(1:n,:));
u=out.model.displacementPressure*pg+out.model.displacementTemperature*(z(4*n+1:end,:)-out.model.referenceTemperature_K);
gap=out.model.gapOperator*u;
a=struct('cycleVolume',sum(volume),'bins',volume,'peakPressure',max(pg,[],'all'), ...
    'peakGap',max(gap,[],'all'),'peakVapor',max(z(2*n+(1:n),:),[],'all'));
end
