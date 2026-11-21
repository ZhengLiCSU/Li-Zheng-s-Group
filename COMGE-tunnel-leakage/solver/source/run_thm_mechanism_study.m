function study=run_thm_mechanism_study(user)
%RUN_THM_MECHANISM_STUDY Independent histories for leakage mechanisms.
% All physical variants start at the same prescribed original T/h forcing,
% and independently equilibrate installation. CP19 is NEVER reused after
% changing ice, conductance, permeability or aperture laws. Defaults request
% full20 histories; maxYears<20 is labelled preliminary, not a Y20 result.
% model may include a fixed soilKsatField realization, shared by all variants.
if nargin<1,user=struct();end
root=fileparts(mfilename('fullpath'));addpath(root);
c=struct('outputFolder',fullfile(root,'leakage_research_20261001','mechanism_histories'), ...
    'maxYears',20,'model',struct(),'resume',false);
names=fieldnames(user);
for k=1:numel(names),assert(isfield(c,names{k}),'MechanismStudy:Config','Unknown option.');c.(names{k})=user.(names{k});end
baseline=joint_contact_thm_model(c.model);
if isfield(c.model,'initial')&&isfield(c.model.initial,'mode')
    assert(any(strcmp(c.model.initial.mode,{'profile','fields'})), ...
        'MechanismStudy:Initial','Independent histories may use a common T/h initial field, not a checkpoint/installation transfer.');
end
variants={'baseline','joint_ice_off','joint_viscosity_off','opening_feedback_off','soil_Ksat_x0p01','joint_Cref_x10'};
overrides={struct(),struct('leakage',struct('jointFreeze',struct('iceEnabled',false))), ...
    struct('leakage',struct('jointFreeze',struct('viscosityEnabled',false))), ...
    struct('thm',struct('openingTransmissionFraction',0)), ...
    struct('Ksat',baseline.Ksat*.01), ...
    struct('leakage',struct('jointConductance',baseline.leakage.jointConductance*10))};
if ~isfolder(c.outputFolder),mkdir(c.outputFolder);end
metrics=cell(size(variants));rows=struct([]);results=cell(size(variants));
for k=1:numel(variants)
    model=combine(c.model,overrides{k});folder=fullfile(c.outputFolder,variants{k});
    fprintf('Independent mechanism history: %s, %d years.\n',variants{k},c.maxYears);
    r=run_TH_segment_year20_thm(struct('outputFolder',folder,'maxYears',c.maxYears, ...
        'model',model,'resume',c.resume));
    assert(r.completedYears==c.maxYears,'MechanismStudy:Incomplete', ...
        'Variant stopped before the requested duration; evidence preserved.');
    metrics{k}=thm_seasonal_leakage_metrics(r.lastYear);results{k}=fullfile(folder,'thm_progress.mat');
    m=metrics{k};rows(k).variant=variants{k};rows(k).completedYears=r.completedYears;
    rows(k).annual_L_m=m.annual_L_m;rows(k).peakDaily_L_m=m.peakDaily_L_m;
    rows(k).firstActiveDay=m.firstActiveDay;rows(k).lastActiveDay=m.lastActiveDay;
    rows(k).signedSameStateSupplyDeficitFraction=m.signedSupplyDeficitFraction;
    save(fullfile(c.outputFolder,'mechanism_study_progress.mat'),'c','metrics','results','rows','-v7');
end
study=struct('complete',true,'isFull20HistoryStudy',c.maxYears==20,'config',c, ...
    'summary',struct2table(rows),'metrics',{metrics},'resultProgressFiles',{results}, ...
    'scope','Independently solved histories; responses include interactions, ablation differences are not additive contributions. Parameters remain illustrative.');
save(fullfile(c.outputFolder,'mechanism_study_report.mat'),'study','-v7');
writetable(study.summary,fullfile(c.outputFolder,'mechanism_summary.csv'));
end
function a=combine(a,b)
f=fieldnames(b);
for k=1:numel(f)
    n=f{k};if isfield(a,n)&&isstruct(a.(n))&&isstruct(b.(n)),a.(n)=combine(a.(n),b.(n));else,a.(n)=b.(n);end
end
end
