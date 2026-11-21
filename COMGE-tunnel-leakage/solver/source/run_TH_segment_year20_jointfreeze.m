function result=run_TH_segment_year20_jointfreeze(user)
%RUN_TH_SEGMENT_YEAR20_JOINTFREEZE Longitudinal-joint TH with reversible icing.
% Default: a complete 20-year eight-segment, one-defect illustrative case.
% The new joint law integrates local microchannel resistance along the real
% lining thickness. It does not introduce cumulative damage or gasket data.
% Soil, climate, water and joint microchannel parameters require calibration.
% Default maximum step: 4 days for annual-volume screening. In the supplied
% example its year-20 volume differs by 0.20% from a 0.5-day final-year
% rerun. Use <=1 day for seasonal transition details; recheck new inputs.
% Example: r=run_TH_segment_year20_jointfreeze;
% Resume:  r=run_TH_segment_year20_jointfreeze(struct('resume',true));
% A smaller maxYears may be supplied for short verification runs; these
% runs deliberately leave the age-20 quantities as NaN.
if nargin<1, user=struct(); end
assert(isstruct(user)&&isscalar(user),'Input must be a scalar struct.');
root=fileparts(mfilename('fullpath')); addpath(root);
if ~isfield(user,'model'), user.model=struct(); end
assert(isstruct(user.model)&&isscalar(user.model),'Model must be a scalar struct.');
model=joint_only_th_model(user.model);
settings=joint_freeze_defaults(struct('enabled',true));
if isfield(model.leakage,'jointFreeze')
    settings=combine(settings,model.leakage.jointFreeze);
end
model.leakage.jointFreeze=joint_freeze_defaults(settings);
if ~isfield(model,'dtMax'), model.dtMax=4*86400; end
assert(isscalar(model.dtMax)&&isfinite(model.dtMax)&& ...
    model.dtMax>0&&model.dtMax<=4*86400, ...
    'The joint-freezing screening maximum step is four days.');
user.model=model;
if ~isfield(user,'outputFolder')
    user.outputFolder=fullfile(root,'segment_year20_jointfreeze_8segment_illustrative');
end
if ~isfield(user,'maxYears'), user.maxYears=20; end
if ~isfield(user,'consecutiveYears'), user.consecutiveYears=21; end
if ~isfield(user,'productionYears'), user.productionYears=1; end
assert(user.maxYears<=20, ...
    'This age-20 wrapper permits at most 20 years; use the general runner otherwise.');
started=tic;
result=run_TH_segment_spinup_v12(user);
result.elapsedWallTime_s=toc(started);
result.jointOnlyModel=true;
result.jointFreezeModel=model.leakage.jointFreeze.enabled;
result.jointParameterStatus=['Eight-segment reference geometry; one assumed ' ...
    'defective joint; reversible microchannel ice blockage and viscosity. ' ...
    'Soil, climate, water and microchannel parameters are uncalibrated.'];
result.year20AnnualLeakage_L_per_m=NaN;
result.year20JointLeakage_L_per_m=NaN;
result.year20BodyLeakage_L_per_m=NaN;
result.year20ByJoint_L_per_m=NaN;
result.year20ByMonth_L_per_m=NaN;
result.year20MonthlyByJoint_L_per_m=NaN;
result.year20DailyLeakage_L_per_day_per_m=NaN;
if result.completedYears==20
    leak=result.lastYear.leakage;
    total=leak.totalVolume_m3_per_m;
    joint=sum(leak.volumeByJoint_m3_per_m);
    body=leak.bodyVolume_m3_per_m;
    assert(result.lastYear.parameters.liningK==0&& ...
        abs(body)<1e-12&&abs(joint-total)<1e-9, ...
        'Joint-only water partition failed.');
    assert(abs(sum(leak.volumeBySegment_m3_per_m)-total)<1e-9, ...
        'Segment rates do not partition total joint leakage.');
    assert(abs(sum(leak.monthlyTotalVolume_m3_per_m)-total)<1e-9, ...
        'Monthly rates do not partition total joint leakage.');
    result.year20AnnualLeakage_L_per_m=1000*total;
    result.year20JointLeakage_L_per_m=1000*joint;
    result.year20BodyLeakage_L_per_m=1000*body;
    result.year20ByJoint_L_per_m=1000*leak.volumeByJoint_m3_per_m;
    result.year20ByMonth_L_per_m=1000*leak.monthlyTotalVolume_m3_per_m;
    result.year20MonthlyByJoint_L_per_m=1000*leak.monthlyVolumeByJoint_m3_per_m;
    result.year20DailyLeakage_L_per_day_per_m=1000*leak.dailyTotalVolume_m3_per_m;
end
% Preserve wrapper metadata after the generic runner's completed-year save.
save(fullfile(user.outputFolder,'TH_segment_spinup_results.mat'),'result','-v7');
end

function a=combine(a,b)
assert(isstruct(b)&&isscalar(b),'Overrides must be scalar structures.');
names=fieldnames(b);
for k=1:numel(names)
    name=names{k};
    if isfield(a,name)&&isstruct(a.(name))&&isstruct(b.(name))
        a.(name)=combine(a.(name),b.(name));
    else
        a.(name)=b.(name);
    end
end
end
