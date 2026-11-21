function report=run_luo2026_constitutive_reproduction(folder)
%RUN_LUO2026_CONSTITUTIVE_REPRODUCTION Public equations, NOT full FE paper.
% Source: Table1, Eqs18,20-26, Fig2; doi:10.1016/j.compgeo.2026.108114.
% Published FE cycle counts are comparison targets, never fitted here.
if nargin<1,folder=fullfile(fileparts(mfilename('fullpath')),'literature_damage_20261002','constitutive_reference');end
if ~isfolder(folder),mkdir(folder);end
s=[.95;.90;.80;.70];[increment,N]=luo2026_fatigue_increment(s);
reportedFE=[8;54;NaN;NaN];
t=table(s,N,ceil(N),increment,reportedFE, ...
    'VariableNames',{'peakStressRatio','equation26Life','materialPointFailureCycle', ...
    'fatigueIncrement','publishedFullFECycle'});
writetable(t,fullfile(folder,'Luo2026_SN_reference.csv'));
material=struct('strength_Pa',9.87e6,'crackingGap_m',1e-7,'debondingGap_m',1e-5);
gap=linspace(0,material.debondingGap_m,401)';sigma=zeros(size(gap));exact=sigma;old=struct();
for k=1:numel(gap)
    [sigma(k),~,old]=luo2026_cohesive_mode1(gap(k),0,old,material);
end
small=gap<=material.crackingGap_m;
exact(small)=material.strength_Pa/material.crackingGap_m*gap(small);
exact(~small)=material.strength_Pa*(material.debondingGap_m-gap(~small))/ ...
    (material.debondingGap_m-material.crackingGap_m);
error=max(abs(sigma-exact))/material.strength_Pa;
writetable(table(gap,sigma,exact),fullfile(folder,'Luo2026_modeI_reference.csv'));
% Eq15, t_l=b^3/12: convert intrinsic transmissivity to head conductance.
% A positive activation tag is used with base mobility zero, so ONLY the
% published ideal slit branch is measured. No calibrated seal is implied.
aperture=logspace(-6,-3,25)';rho=1000;gravity=9.81;length=.45;
options=struct('openingTransmissionFraction',1,'baseMobilityFactor',0, ...
    'pathLength_m',length,'waterDensity_kg_m3',rho,'gravity_m_s2',gravity, ...
    'referenceTemperature_C',20,'freezeEnabled',false,'viscosityEnabled',false);
C=joint_opening_hydraulics(ones(size(aperture))*1e-10,aperture,[],options);
exactC=rho*gravity*aperture.^3/(12*water_viscosity_iapws(20)*length);
hydraulicError=max(abs(C-exactC)./exactC);
writetable(table(aperture,C,exactC),fullfile(folder,'Luo2026_cubic_hydraulic_reference.csv'));
% Eq25 differs from a bare sum(1/N). These are law-only predictions.
paths={.95*ones(100,1),.9*ones(100,1),[.9*ones(10,1);.95*ones(90,1)], ...
    [.95*ones(5,1);.9*ones(95,1)]};names={'constant95','constant90','lowThenHigh','highThenLow'};
rows=[];summary=zeros(4,2);
for j=1:4
    [inc,~]=luo2026_fatigue_increment(paths{j});D=cumsum(inc);
    failure=find(1-D<=paths{j}+1e-14,1);summary(j,:)=[j,failure];
    rows=[rows;repmat(j,100,1),(1:100)',paths{j},D,1-D]; %#ok<AGROW>
end
writetable(array2table(rows,'VariableNames',{'caseID','cycle','peakRatio','Df','strengthRatio'}), ...
    fullfile(folder,'Luo2026_variable_amplitude.csv'));
report=struct('passed',error<1e-12&&hydraulicError<1e-12,'relativeTractionError',error, ...
    'relativeCubicHydraulicError',hydraulicError, ...
    'SN',t,'variableAmplitudeFailureCycles',summary,'caseNames',{names}, ...
    'source','https://doi.org/10.1016/j.compgeo.2026.108114', ...
    'scope','Equation/material-point reproduction only; no Brazilian mesh, fracture network or experimental reproduction.', ...
    'discrepancy','Table1 S-N law gives ceil(life)=7/66 for 95%/90%, while the paper reports full-FE failure at 8/54. FE local stress, viscosity and discretization are not reproduced. Fig10 text also uses both 14 and 15 for CaseA; no forced agreement.');
save(fullfile(folder,'Luo2026_reference_report.mat'),'report','material','-v7');
assert(report.passed,'Luo2026:Reproduction','Mode-I constitutive identity failed.');
end
