function report = run_premelting_verification(folder)
%RUN_PREMELTING_VERIFICATION Reproducible constitutive + inventory campaign.
% Keep this evidence separate from the validated full THM 20-year branch.
arguments
    folder (1,:) char
end
assert(~isfolder(folder),'Premelting:OutputExists','Use a fresh folder; preserve prior evidence.');
mkdir(folder); mkdir(fullfile(folder,'source'));
files=dir('*.m');
for k=1:numel(files), copyfile(files(k).name,fullfile(folder,'source',files(k).name)); end
tests=runtests({'TestJointThermodynamics.m','TestNonisothermalInterface.m','TestThreePhaseInterface.m', ...
    'TestCohesiveThreePhase.m','TestPremelting.m'});
save(fullfile(folder,'tests.mat'),'tests','-v7');
assert(all([tests.Passed]),'Premelting:Tests','Verification tests did not all pass.');
p=joint_threephase_parameters(); f=premelting_film_parameters();
T=linspace(273.05,268.15,121)'; pl=2e5*ones(size(T)); aperture=5e-6*ones(size(T));
[s,d]=premelting_slit_state(pl,T,aperture,p,f);
identity=d.freeEnergy+s.entropy_Jm2K.*[0,1,0]+s.equilibrium.icePressure_Pa.*[0,0,1]- ...
    s.equilibrium.chemicalPotential_Jkg.*d.mass;
relativeIdentity=max(abs(identity)./max(1,abs(d.freeEnergy)),[],'all');
assert(relativeIdentity<1e-10,'Premelting:EnergyIdentity','Free-energy differential identity failed.');
coexistence=zeros(4,3);
for j=1:4
    material=premelting_film_parameters(j);
    root=fzero(@(h)coexistencePressure(h,material),[material.minimumFilm_m,material.maximumFilm_m],optimset('TolX',1e-18));
    coexistence(j,:)=[j,root,material.tabulatedCoexistenceFilm_m];
end
% Published table thicknesses have only 2--3 significant digits; rounded
% fitted coefficients need not reproduce them at machine precision.
coexistenceError=max(abs(coexistence(:,2)-coexistence(:,3))./coexistence(:,3));
assert(coexistenceError<0.02,'Premelting:Literature','Published film-thickness reproduction exceeds 2 percent.');
save(fullfile(folder,'constitutive.mat'),'s','d','p','f','coexistence','identity','-v7');

% Mechanical OPEN-RESERVOIR equilibrium limit: no calendar time is assigned.
% Local pl is imposed everywhere, so changing inventory assumes replenishment.
c=joint_cohesive_parameters(struct('strength_Pa',3e5));
m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
n=numel(m.weight); Ts=m.referenceTemperature_K*ones(numel(m.capacity),1);
loadingT=linspace(273.05,270.15,61); pathT=[loadingT,fliplr(loadingT(1:end-1))];
history=zeros(n,1); path=cell(numel(pathT),1);
for k=1:numel(pathT)
    path{k}=premelting_ring_state(2e5*ones(n,1),pathT(k)*ones(n,1),Ts,m,p,f,history);
    history=path{k}.cohesive.kappa_m;
end
before=cohesive_isothermal_leak_probe(m,p,zeros(n,1));
after=cohesive_isothermal_leak_probe(m,p,history);
save(fullfile(folder,'open_equilibrium_limit.mat'),'m','p','f','pathT','path','before','after','-v7');

coarse=run_premelting_inventory(1,false); save(fullfile(folder,'inventory_day.mat'),'coarse','-v7');
fine=run_premelting_inventory(0.5,false); save(fullfile(folder,'inventory_halfday.mat'),'fine','-v7');
quarter=run_premelting_inventory(0.25,false); save(fullfile(folder,'inventory_quarterday.mat'),'quarter','-v7');
sealed=run_premelting_inventory(1,true); save(fullfile(folder,'inventory_sealed.mat'),'sealed','-v7');
report=struct('testsPassed',sum([tests.Passed]),'coexistenceRelativeError',coexistenceError, ...
    'freeEnergyIdentityRelativeError',relativeIdentity, ...
    'scope','Constitutive, Q4 equilibrium and prescribed-temperature frozen inventory verification only.', ...
    'fullTransientTHMValidated',false,'new20YearRunValidated',false);
save(fullfile(folder,'verification_matlab.mat'),'report','-v7');
end

function r=coexistencePressure(h,f)
v=premelting_film_potential(h,f); r=v.pressure_Pa;
end
