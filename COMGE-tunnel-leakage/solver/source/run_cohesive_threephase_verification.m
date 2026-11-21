function files = run_cohesive_threephase_verification(folder)
%RUN_COHESIVE_THREEPHASE_VERIFICATION Mechanism controls, NOT tunnel years.
% A weak BONDED interface example; not an EPDM material recommendation.
% Fracture work is retained as surface energy, not an empirical heat source.
arguments
    folder (1,:) char
end
assert(~isfolder(folder),'Cohesive:OutputExists','Use a fresh directory to preserve evidence.');
mkdir(folder); source=fullfile(folder,'source'); mkdir(source);
here=fileparts(mfilename('fullpath')); list=dir(fullfile(here,'*.m'));
for k=1:numel(list), copyfile(fullfile(list(k).folder,list(k).name),source); end
p=joint_threephase_parameters();
c=joint_cohesive_parameters(struct('strength_Pa',3e5));
S=make_thm_verification_mesh(16);
isolated=build_cohesive_threephase_ring(S,3,p,c,struct('thermalExpansion_Kinv',0));
normal=build_cohesive_threephase_ring(S,3,p,c);
elastic=c; elastic.enabled=false;
control=build_cohesive_threephase_ring(S,3,p,elastic,struct('thermalExpansion_Kinv',0));
names={'phase_expansion_damage_reference','phase_expansion_damage_fine', ...
    'phase_expansion_damage_quarter','phase_expansion_damage_eighth','phase_expansion_elastic_control','warm_damage_control','thermal_damage_reference','thermal_damage_fine'};
models={isolated,isolated,isolated,isolated,control,isolated,normal,normal};
base=struct('cycles',2,'temperatureAmplitude_K',20);
files=cell(size(names));
for k=1:numel(names)
    o=base;
    if (k>=2 && k<=4) || k==8
        factor=2^(k-1);
        if k==8, factor=2; end
        o.dtDays=0.125/factor; o.maximumPhaseIncrement=0.025/factor;
        o.maximumPressureIncrement_Pa=1e4/factor; o.maximumLogPressureIncrement=0.25/factor;
        o.maximumDamageIncrement=0.025/factor;
    elseif k==6
        o.temperatureAmplitude_K=0;
    end
    files{k}=fullfile(folder,[names{k} '.mat']); o.saveFile=files{k};
    solve_cohesive_threephase_ring(models{k},p,o);
end
fprintf('Completed mechanism runs. Independently audit using audit_cohesive_threephase.py.\n');
end
