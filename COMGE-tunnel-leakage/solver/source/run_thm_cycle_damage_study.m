function study=run_thm_cycle_damage_study(user)
%RUN_THM_CYCLE_DAMAGE_STUDY Independent continuous local-cycle histories.
% Default strict 20 years; maxYears<20 is explicitly only preliminary.
% doseRates compare a zero-damage control and uncalibrated coarsening rates.
% Does not borrow CP19, count calendar years, or soften gasket mechanics.
if nargin<1,user=struct();end
root=fileparts(mfilename('fullpath'));
c=struct('outputFolder',fullfile(root,'literature_damage_20261002','strict_cycle20'), ...
    'maxYears',20,'doseRates',[0 .05 .2],'cycleSettings',struct(), ...
    'model',struct(),'resume',false,'useFastNumerics',false);
for key=fieldnames(user)'
    assert(isfield(c,key{1}),'JointCycleStudy:Config','Unknown option: %s.',key{1});
    c.(key{1})=user.(key{1});
end
assert(isnumeric(c.doseRates)&&isreal(c.doseRates)&&isvector(c.doseRates) ...
    &&~isempty(c.doseRates)&&all(isfinite(c.doseRates(:))) ...
    &&all(c.doseRates(:)>=0)&&numel(unique(c.doseRates))==numel(c.doseRates), ...
    'JointCycleStudy:Config','Supply distinct nonnegative dose rates.');
assert(c.maxYears>=1&&c.maxYears<=20&&fix(c.maxYears)==c.maxYears ...
    &&isscalar(c.maxYears),'JointCycleStudy:Config','maxYears must be 1..20.');
folder=char(java.io.File(c.outputFolder).getCanonicalPath());c.outputFolder=folder;
if ~isfolder(folder),mkdir(folder);end
configFile=fullfile(folder,'cycle_study_configuration.mat');
if c.resume
    old=load(configFile,'configuration');
    assert(isequaln(rmfield(c,'resume'),rmfield(old.configuration,'resume')), ...
        'JointCycleStudy:Resume','Resume requires the identical study configuration.');
else
    assert(~isfile(configFile),'JointCycleStudy:Exists','Use resume=true or a fresh folder.');
    configuration=c;save(configFile,'configuration','-v7');
    mkdir(fullfile(folder,'source_at_start'));
    src=dir(fullfile(root,'*.m'));
    for k=1:numel(src),copyfile(fullfile(root,src(k).name),fullfile(folder,'source_at_start',src(k).name));end
end
critical={'TH_2D_segment_v16.m','joint_freeze_thaw_memory.m','joint_cycle_channel_law.m', ...
    'joint_cycle_damage_defaults.m','joint_cycle_severity.m','validate_joint_cycle_memory.m','joint_channel_inventory.m', ...
    'audit_joint_channel_geometry.m', ...
    'joint_freeze_conductance.m','thm_constitutive_signature.m','frozen_soil_constitutive.m', ...
    'thm_bulk_constitutive.m','joint_opening_p1_hydraulics.m','solve_segment_contact_fem.m', ...
    'run_TH_segment_year20_thm.m','run_TH_segment_year20_fast_thm.m', ...
    'soil_cvfe_flux.m','soil_cvfe_heat.m','soil_poroelastic_fem.m','evaluate_joint_drain.m', ...
    'joint_microchannel_damage.m','joint_opening_profile_hydraulics.m'};
for k=1:numel(critical)
    assert(strcmp(fileread(fullfile(root,critical{k})), ...
        fileread(fullfile(folder,'source_at_start',critical{k}))), ...
        'JointCycleStudy:Source','Physical source changed after the study started: %s.',critical{k});
end
study=struct('configuration',c,'isFull20HistoryStudy',c.maxYears==20,'cases',{{}}, ...
    'summary',table(),'scope','Illustrative constant-volume channel coarsening; numerical verification, no experimental calibration.');
rows=zeros(numel(c.doseRates),8);
for j=1:numel(c.doseRates)
    settings=joint_cycle_damage_defaults(c.cycleSettings);
    settings.enabled=true;settings.doseRate=c.doseRates(j);
    model=c.model;
    if ~isfield(model,'leakage'),model.leakage=struct();end
    model.leakage.jointCycleDamage=settings;
    caseFolder=fullfile(folder,sprintf('rate_%02d',j));
    resumeCase=c.resume&&isfile(fullfile(caseFolder,'thm_progress.mat'));
    cfg=struct('maxYears',c.maxYears,'outputFolder',caseFolder,'model',model,'resume',resumeCase);
    if c.useFastNumerics,result=run_TH_segment_year20_fast_thm(cfg);
    else,result=run_TH_segment_year20_thm(cfg);end
    [audit,~]=validate_thm_20year(caseFolder,struct('makePlots',false));
    yearly=audit.checks(audit.checks.year>0,:);
    assert(all(yearly.passed),'JointCycleStudy:Audit','An original numerical/physical gate failed.');
    out=result.lastYear;memoryAudit=audit_joint_cycle_damage(out);
    assert(memoryAudit.passed,'JointCycleStudy:Audit','Accepted local-cycle history failed.');
    rows(j,:)=[c.doseRates(j),c.maxYears,1000*out.leakage.totalVolume_m3_per_m, ...
        out.jointCycleDamage.appliedSeverity(7),mean(out.jointCycleDamage.finalMemory.cycleCount(7,:)), ...
        out.solverWork.acceptedSteps,out.solverWork.rejectedSteps, ...
        max(abs(out.history(:,8)))];
    study.cases{j}=struct('folder',caseFolder,'audit',audit,'memoryAudit',memoryAudit);
    h=out.jointCycleDamage.history;
    writetable(array2table([out.history(:,1),out.history(:,10)*86400*1000,h(:,8),h(:,24),h(:,40),h(:,48)], ...
        'VariableNames',{'day','totalLeak_L_per_day_m','J7AppliedSeverity', ...
        'J7MeanCycles','J7ConnectedVolume_m2','J7Cref_m_s'}),fullfile(caseFolder,'last_year_cycle_daily_rate.csv'));
    study.summary=array2table(rows(1:j,:), ...
        'VariableNames',{'doseRate','completedYears','lastYearLeak_L_m','J7AppliedSeverity', ...
        'J7MeanLocalCycles','lastYearAcceptedSteps','lastYearRejectedSteps','maxWaterDefect_m2'});
    save(fullfile(folder,'cycle_study_report.mat'),'study','-v7');
    writetable(study.summary,fullfile(folder,'cycle_study_summary.csv'));
end
end
