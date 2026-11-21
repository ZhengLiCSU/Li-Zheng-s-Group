function study=run_thm_joint_damage_study(user)
%RUN_THM_JOINT_DAMAGE_STUDY Independent THM histories at fixed joint damage.
% Starts each geometry case at the prescribed original T/h, independently
% equilibrates installation, and retains all requested consecutive years.
% Different damage cases never borrow CP19 or reset inventories mid-run.
% This is a damage-state sensitivity study, NOT a cycle-to-damage law.
% Default strict mesh, 20 years. Shorter and fast cases are labelled as such.
arguments
    user (1,1) struct = struct()
end
root=fileparts(mfilename('fullpath'));
c=struct('outputFolder',fullfile(root,'joint_damage_20261002','strict_full20'), ...
    'levels',[0 .5 1],'maxYears',20,'model',struct(),'resume',false, ...
    'radiusGrowthAtUnitDamage',.5,'channelDensityGrowthAtUnitDamage',0, ...
    'damagePattern',[],'useFastNumerics',false);
keys=fieldnames(user);
for k=1:numel(keys)
    assert(isfield(c,keys{k}),'JointDamageStudy:Config','Unknown option: %s.',keys{k});
    c.(keys{k})=user.(keys{k});
end
assert(isnumeric(c.levels)&&isreal(c.levels)&&isvector(c.levels)&&all(isfinite(c.levels)) ...
    &&all(c.levels>=0&c.levels<=1)&&all(diff(c.levels)>0), ...
    'JointDamageStudy:Levels','Use strictly increasing damage levels in [0,1].');
assert(isscalar(c.useFastNumerics)&&ismember(c.useFastNumerics,[0 1]), ...
    'JointDamageStudy:Config','useFastNumerics must be logical 0/1.');
baseline=joint_contact_thm_model(c.model);
if isfield(baseline,'initial')&&isfield(baseline.initial,'mode')
    assert(any(strcmp(baseline.initial.mode,{'profile','fields'})), ...
        'JointDamageStudy:Initial','Changed damage cases require independent original histories.');
end
nj=numel(baseline.leakage.jointAnglesRad);
baseC=baseline.leakage.jointConductance(:);
if isscalar(baseC),baseC=repmat(baseC,nj,1);end
if isempty(c.damagePattern),c.damagePattern=double(baseC>0);end
assert(isnumeric(c.damagePattern)&&isreal(c.damagePattern)&&numel(c.damagePattern)==nj ...
    &&all(isfinite(c.damagePattern(:)))&&all(c.damagePattern(:)>=0&c.damagePattern(:)<=1), ...
    'JointDamageStudy:Pattern','Supply one normalized damage susceptibility per joint.');
folderObject=java.io.File(c.outputFolder);c.outputFolder=char(folderObject.getCanonicalPath());
if ~isfolder(c.outputFolder),mkdir(c.outputFolder);end
configurationFile=fullfile(c.outputFolder,'damage_study_configuration.mat');
if isfile(configurationFile)
    saved=load(configurationFile,'c');
    assert(c.resume&&isequaln(rmfield(saved.c,'resume'),rmfield(c,'resume')), ...
        'JointDamageStudy:Config','Existing study requires identical physical cases and numerical configuration.');
else
    save(configurationFile,'c','-v7');
end
sourceFolder=fullfile(c.outputFolder,'source_at_start');
if ~isfolder(sourceFolder)
    mkdir(sourceFolder);files=dir(fullfile(root,'*.m'));
    for k=1:numel(files),copyfile(fullfile(root,files(k).name),fullfile(sourceFolder,files(k).name));end
else
    assert(c.resume,'JointDamageStudy:Existing','Existing evidence requires resume=true or a new folder.');
    critical={'TH_2D_segment_v16.m','joint_microchannel_damage.m','thm_constitutive_signature.m', ...
        'frozen_soil_constitutive.m','thm_bulk_constitutive.m','joint_channel_inventory.m', ...
        'joint_freeze_conductance.m','evaluate_joint_drain.m','solve_segment_contact_fem.m', ...
        'joint_opening_p1_hydraulics.m','soil_cvfe_flux.m','soil_cvfe_heat.m', ...
        'soil_poroelastic_fem.m','run_TH_segment_year20_thm.m','run_TH_segment_year20_fast_thm.m'};
    for k=1:numel(critical)
        assert(strcmp(fileread(fullfile(sourceFolder,critical{k})),fileread(fullfile(root,critical{k}))), ...
            'JointDamageStudy:Source','Execution source differs from the frozen source: %s.',critical{k});
    end
end
metrics=cell(numel(c.levels),1);rows=struct([]);
for k=1:numel(c.levels)
    model=c.model;
    if ~isfield(model,'leakage'),model.leakage=struct();end
    model.leakage.jointDamage=struct('enabled',true, ...
        'severity',c.levels(k)*c.damagePattern(:), ...
        'radiusGrowthAtUnitDamage',c.radiusGrowthAtUnitDamage, ...
        'channelDensityGrowthAtUnitDamage',c.channelDensityGrowthAtUnitDamage);
    caseFolder=fullfile(c.outputFolder,sprintf('damage_%02d',k));
    caseResume=c.resume&&isfile(fullfile(caseFolder,'thm_progress.mat'));
    cfg=struct('outputFolder',caseFolder,'model',model,'maxYears',c.maxYears,'resume',caseResume);
    fprintf('Fixed residual joint damage %.4g; independent %d-year history.\n',c.levels(k),c.maxYears);
    if c.useFastNumerics
        r=run_TH_segment_year20_fast_thm(cfg);
    else
        r=run_TH_segment_year20_thm(cfg);
    end
    assert(r.completedYears==c.maxYears,'JointDamageStudy:Incomplete','Incomplete evidence retained.');
    audit=validate_thm_20year(caseFolder,struct('makePlots',false,'saveReport',false));
    assert(all(audit.checks.passed(audit.checks.year>0)), ...
        'JointDamageStudy:Audit','A completed-year physical/continuity audit failed.');
    save(fullfile(caseFolder,'damage_completed_years_audit.mat'),'audit','-v7');
    m=thm_seasonal_leakage_metrics(r.lastYear);metrics{k}=m;
    writetable(m.daily,fullfile(c.outputFolder,sprintf('damage_%02d_daily.csv',k)));
    writetable(table((1:12)',m.monthly_L_m,'VariableNames',{'month','leakage_L_m'}), ...
        fullfile(c.outputFolder,sprintf('damage_%02d_monthly.csv',k)));
    g=r.lastYear.jointDamage;vol=r.lastYear.state.jointInventory.waterVolume_m2;
    rows(k).damageLevel=c.levels(k);rows(k).completedYears=r.completedYears;
    rows(k).annual_L_m=m.annual_L_m;rows(k).peakDaily_L_m=m.peakDaily_L_m;
    rows(k).firstActiveDay=m.firstActiveDay;rows(k).lastActiveDay=m.lastActiveDay;
    rows(k).radiusMultiplier_activeJoint=max(g.radiusMultiplier(baseC>0));
    rows(k).referenceConductanceRatio_activeJoint=max(g.referenceConductanceMultiplier(baseC>0));
    rows(k).connectedWaterVolume_m2=sum(vol);
    rows(k).signedSupplyDeficitFraction=m.signedSupplyDeficitFraction;
    rows(k).bodyLeak_m3_m=r.lastYear.leakage.bodyVolume_m3_per_m;
    rows(k).totalAcceptedSteps=sum(r.annual(:,9));
    rows(k).lastYearRejectedSteps=r.lastYear.solverWork.rejectedSteps;
    rows(k).caseFolder=caseFolder;
    save(fullfile(c.outputFolder,'damage_study_progress.mat'),'c','rows','metrics','-v7');
end
study=struct('complete',true,'isFull20HistoryStudy',c.maxYears==20,'config',c, ...
    'summary',struct2table(rows),'metrics',{metrics}, ...
    'annualLeakageIsIncreasing',all(diff([rows.annual_L_m])>0), ...
    'scope','Fixed damaged joint-interface geometries, independent histories; not freeze-thaw accumulation, EPDM constitutive damage, or calibrated field prediction.');
save(fullfile(c.outputFolder,'joint_damage_study_report.mat'),'study','-v7');
writetable(study.summary,fullfile(c.outputFolder,'damage_summary.csv'));
end
