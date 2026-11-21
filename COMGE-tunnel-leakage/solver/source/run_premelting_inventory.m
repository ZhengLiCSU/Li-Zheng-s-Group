function out = run_premelting_inventory(dtDays,sealed)
%RUN_PREMELTING_INVENTORY Temperature-controlled FROZEN-ring benchmark.
% Two 20-day cooling/rewarming excursions of 0.1 K; ice remains throughout.
% This is NOT two complete freezing/thawing cycles or a seasonal heat solve.
arguments
    dtDays (1,1) double {mustBePositive,mustBeFinite} = 1
    sealed (1,1) logical = false
end
p=joint_threephase_parameters(); f=premelting_film_parameters();
c=joint_cohesive_parameters(struct('enabled',false));
m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
n=numel(m.weight); ns=numel(m.capacity); pl=2e5*ones(n,1); T0=272.15;
Ts=m.referenceTemperature_K*ones(ns,1);
s=premelting_ring_state(pl,T0*ones(n,1),Ts,m,p,f,zeros(n,1));
bc=thermo_interface_boundary(m,279.15,279.15,2e5,sealed);
timeDays=unique([0:dtDays:40,20,40])'; nt=numel(timeDays);
states=cell(nt,1); states{1}=s; audits=cell(nt-1,1);
history=zeros(n,nt); history(:,1)=s.cohesive.kappa_m;
for k=2:nt
    T=(T0-0.1*sin(pi*timeDays(k)/20)^2)*ones(n,1);
    [s,a]=premelting_prescribed_temperature_step(s,T,(timeDays(k)-timeDays(k-1))*86400,m,bc,p,f);
    states{k}=s; audits{k-1}=a; history(:,k)=s.cohesive.kappa_m;
end
out=struct('engine','PremeltingInventory_v1_prescribed_temperature','timeDays',timeDays, ...
    'states',{states},'audit',{audits},'history',history,'model',m,'parameters',p,'filmParameters',f, ...
    'boundary',bc,'sealed',sealed,'completed',true, ...
    'scope','Mass-conservative Q4 wet-ice solve with prescribed temperatures; no heat PDE, ice nucleation or melting.');
end
