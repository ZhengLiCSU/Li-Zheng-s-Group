function study=run_thm_seasonal_convergence(user)
%RUN_THM_SEASONAL_CONVERGENCE Conditional Y20 path/time refinement.
% Complete source history is immutable. Each refinement starts from the
% SAME full CP19 and installed U/T. Only time/path/SCV quadrature settings
% may change, with explicit new-run reintegration; no saved history reused.
% Default source is strict 20y. levels rows: [dtMax_days radialSamples].
% Default refines radial resolution at fixed dt, then dt at fixed samples.
% Study compares against the last level. All differences are reported;
% no annual tolerance is used to hide peak/date/small-month discrepancies.
if nargin<1,user=struct();end
root=fileparts(mfilename('fullpath'));addpath(root);
c=struct('sourceFolder',fullfile(root,'thm20_20260930_analytic_p1'), ...
    'outputFolder',fullfile(root,'leakage_research_20261001','seasonal_convergence'), ...
    'levels',[],'activityThresholdFraction',.01);
c=mergeConfig(c,user);if ~isfolder(c.outputFolder),mkdir(c.outputFolder);end
d=load(fullfile(c.sourceFolder,'thm_progress.mat'),'progress');p=d.progress;
assert(p.completedYears==20&&numel(p.coreCheckpoints)==20, ...
    'SeasonalStudy:Source','A registered complete20 history is required.');
base=p.lastYear;cpFile=fullfile(c.sourceFolder,'year_19_checkpoint.mat');
cp=load(cpFile,'checkpoint');cp=cp.checkpoint;
assert(isequal(base.initialX,cp.X)&&isequal(cp.mechanicalReference.U,base.mechanics.referenceU) ...
    &&isequal(cp.mechanicalReference.T,base.mechanics.referenceTemperature), ...
    'SeasonalStudy:Source','Saved target-year X and installed reference disagree.');
if isempty(c.levels)
    dt=base.parameters.dtMax/86400;n=base.parameters.leakage.jointFreeze.quadraturePoints;
    c.levels=[dt n;dt 32;dt 64;dt/2 64;dt/4 64];
end
assert(isnumeric(c.levels)&&size(c.levels,2)==2&&size(c.levels,1)>=2 ...
    &&all(isfinite(c.levels(:)))&&all(c.levels(:,1)>0) ...
    &&all(c.levels(:,2)>=1&c.levels(:,2)==fix(c.levels(:,2))), ...
    'SeasonalStudy:Levels','Supply at least two [dtMax_days radialSamples] rows.');
metrics=cell(size(c.levels,1),1);files=cell(size(metrics));row=struct([]);
for k=1:numel(metrics)
    prefix=fullfile(c.outputFolder,sprintf('level_%02d',k));file=[prefix '_results.mat'];
    o=base.parameters;o.mesh=base.mesh;o.durationDays=365;o.dtMax=c.levels(k,1)*86400;
    o.dtInitial=min(o.dtInitial,o.dtMax);o.leakage.jointFreeze.quadraturePoints=c.levels(k,2);
    % Scale the original calendar caps with dt; otherwise an advertised time
    % refinement could leave already-smaller transition-season caps unchanged.
    if ~isempty(o.solver.dtMaxByMonth_days)
        o.solver.dtMaxByMonth_days=min(o.solver.dtMaxByMonth_days, ...
            o.solver.dtMaxByMonth_days*c.levels(k,1)/(base.parameters.dtMax/86400));
    end
    o.initial=struct('mode','checkpoint','checkpointFile',cpFile, ...
        'resumeHistory',false,'allowNumericalReintegration',true);
    o.outputPrefix=prefix;o.plot=false;o.verbose=false;o.saveOutput=true;
    if isfile(file)
        saved=load(file,'out','studyOptions');
        assert(isfield(saved,'studyOptions')&&isequaln(saved.studyOptions,o), ...
            'SeasonalStudy:Existing','Existing file has a different source/configuration; preserved.');
        out=saved.out;
    else
        fprintf('Seasonal conditional Y20 level %d/%d: dt<=%.4g d, radial samples=%d.\n', ...
            k,numel(metrics),c.levels(k,1),c.levels(k,2));
        out=TH_2D_segment_v16(o);studyOptions=o;save(file,'studyOptions','-append');
    end
    assert(isequal(out.initialX,cp.X)&&isequal(out.mechanics.referenceU,cp.mechanicalReference.U) ...
        &&isequal(out.mechanics.referenceTemperature,cp.mechanicalReference.T) ...
        &&isequaln(out.mesh,base.mesh)&&out.leakage.bodyVolume_m3_per_m==0 ...
        &&thm_numerical_reintegration_compatible(thm_constitutive_signature(out.parameters),cp.signature), ...
        'SeasonalStudy:Invariant','Refinement changed starting state/reference/mesh/physical law.');
    metrics{k}=thm_seasonal_leakage_metrics(out,c.activityThresholdFraction);files{k}=file;
    m=metrics{k};row(k).level=k;row(k).dtMax_days=c.levels(k,1);row(k).radialSamples=c.levels(k,2);
    row(k).annual_L_m=m.annual_L_m;row(k).peakDaily_L_m=m.peakDaily_L_m;
    row(k).peakDay=m.peakDailyDay;row(k).firstActiveDay=m.firstActiveDay;
    row(k).lastActiveDay=m.lastActiveDay;row(k).acceptedSteps=m.acceptedSteps;row(k).rejectedSteps=m.rejectedSteps;
    save(fullfile(c.outputFolder,'seasonal_study_progress.mat'),'c','metrics','files','row','-v7');
end
comparison=cellfun(@(m)compare_thm_seasonal_metrics(m,metrics{end}),metrics,'UniformOutput',false);
study=struct('complete',true,'config',c,'metrics',{metrics},'comparisonToLastLevel',{comparison}, ...
    'summary',struct2table(row),'resultFiles',{files}, ...
    'scope','Conditional Y20 spatial path/time refinement from identical CP19; no full20 history or soil mesh refinement, and no experimental validation.');
save(fullfile(c.outputFolder,'seasonal_convergence_report.mat'),'study','-v7');
writetable(study.summary,fullfile(c.outputFolder,'seasonal_convergence_summary.csv'));
for k=1:numel(metrics)
    writetable(metrics{k}.daily,fullfile(c.outputFolder,sprintf('level_%02d_daily.csv',k)));
end
end
function c=mergeConfig(c,user)
assert(isstruct(user)&&isscalar(user),'SeasonalStudy:Config','Config must be a scalar struct.');
f=fieldnames(user);
for k=1:numel(f)
    assert(isfield(c,f{k}),'SeasonalStudy:Config','Unknown field %s.',f{k});c.(f{k})=user.(f{k});
end
end
