function study=run_thm_literature_parameter_study(user)
%RUN_THM_LITERATURE_PARAMETER_STUDY Continuous independent literature-soil cases.
% Defaults to 20 strict years, using the current local-cycle hydraulic law.
% Soil data are traceable; cycle/joint/boundary/mechanical values remain
% separate scenario inputs, not extracted from rock or concrete fatigue.
if nargin<1,user=struct();end
root=fileparts(mfilename('fullpath'));
c=struct('parameterPack','kanagawa_th','maxYears',20,'doseRates',[0 .05 .2], ...
    'cycleSettings',struct(),'model',struct(),'resume',false,'useFastNumerics',false, ...
    'outputFolder',fullfile(root,'literature_parameters_20261002','kanagawa_strict20'));
assert(isstruct(user)&&isscalar(user),'THMLiteratureStudy:Config','Supply scalar struct.');
for key=fieldnames(user)'
    assert(isfield(c,key{1}),'THMLiteratureStudy:Config','Unknown option: %s.',key{1});
    c.(key{1})=user.(key{1});
end
[model,soilEvidence]=thm_literature_parameter_pack(c.parameterPack);
model=combine(model,c.model);
% Reconcile geothermal gradient with changed unfrozen saturated bulk k.
% Boundary flux/air temperature are still the original assumed tunnel case.
ks=2.3;if isfield(model,'ks'),ks=model.ks;end
kw=.57;if isfield(model,'kw'),kw=model.kw;end
geo=.06;if isfield(model,'geothermalFlux'),geo=model.geothermalFlux;end
gradient=geo/((1-model.porosity)*ks+model.porosity*kw);
if ~isfield(model,'initial'),model.initial=struct();end
if ~isfield(model.initial,'gradient'),model.initial.gradient=gradient;end
if ~isfield(model,'groundwaterInflowGradient'),model.groundwaterInflowGradient=gradient;end
soilEvidence.resolvedOverrides=model;
soilEvidence.gradientAssumptions=struct('geothermalFlux_W_m2',geo, ...
    'solidConductivity_W_mK',ks,'waterConductivity_W_mK',kw, ...
    'derivedGradient_C_m',gradient,'source','assumed tunnel case, not soil-column measurements');
soilEvidence.cycleInputs=struct('doseRates',c.doseRates,'settings',c.cycleSettings, ...
    'source','explicit scenario or caller-supplied calibration, not inferred by soil pack');
effective=zeros(height(soilEvidence.ledger),1);
for k=1:height(soilEvidence.ledger)
    key=char(soilEvidence.ledger.parameter(k));
    if strcmp(key,'rhoS*cs'),effective(k)=model.rhoS*model.cs;
    else,effective(k)=model.(key);end
end
soilEvidence.ledger.effectiveValue=effective;
soilEvidence.ledger.overridden=effective~=soilEvidence.ledger.modelValue;
folder=char(c.outputFolder);
if ~java.io.File(folder).isAbsolute(),folder=fullfile(pwd,folder);end
folder=char(java.io.File(folder).getCanonicalPath());
if ~isfolder(folder),mkdir(folder);end
file=fullfile(folder,'parameter_pack_evidence.mat');
if isfile(file)
    old=load(file,'soilEvidence');
    assert(c.resume&&isequaln(old.soilEvidence,soilEvidence), ...
        'THMLiteratureStudy:Resume','Existing case requires identical provenance/parameters and resume=true.');
else
    save(file,'soilEvidence','-v7');
    writetable(soilEvidence.ledger,fullfile(folder,'soil_parameter_ledger.csv'));
end
cfg=struct('outputFolder',folder,'maxYears',c.maxYears,'doseRates',c.doseRates, ...
    'cycleSettings',c.cycleSettings,'model',model,'resume',c.resume, ...
    'useFastNumerics',c.useFastNumerics);
study=run_thm_cycle_damage_study(cfg);study.soilEvidence=soilEvidence;
save(fullfile(folder,'literature_parameter_study.mat'),'study','-v7');
end
function a=combine(a,b)
assert(isstruct(b)&&isscalar(b),'THMLiteratureStudy:Model','Overrides must be scalar struct.');
for key=fieldnames(b)'
    name=key{1};
    if isfield(a,name)&&isstruct(a.(name))&&isstruct(b.(name))
        a.(name)=combine(a.(name),b.(name));
    else
        a.(name)=b.(name);
    end
end
end
