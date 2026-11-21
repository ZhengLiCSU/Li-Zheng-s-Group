function result=run_TH_segment_year20_regime(user)
%RUN_TH_SEGMENT_YEAR20_REGIME Joint-freeze TH with water-supply diagnostics.
% Default: complete 20-year history of the illustrative eight-segment model,
% starting from its initial temperature profile and hydrostatic water head.
% Soil viscosity and finite microchannel water/energy storage are enabled.
% The model does not include cumulative gasket damage or a site calibration.
% Each parameter case needs its own output folder and full thermal/water
% history. A changed parameter case must not resume an old model checkpoint.
%
% r=run_TH_segment_year20_regime;
% r=run_TH_segment_year20_regime(struct('resume',true));
% Partial MODEL overrides are combined recursively with the joint baseline.
% Set maxYears<20 only for short checks; year-20 response fields remain NaN.
% Maximum step is four days for screening and must be independently checked
% for new soil, joint and seasonal water-head scenarios.
if nargin<1, user=struct(); end
assert(isstruct(user)&&isscalar(user),'Input must be a scalar struct.');
root=fileparts(mfilename('fullpath')); addpath(root);
if ~isfield(user,'model'), user.model=struct(); end
assert(isstruct(user.model)&&isscalar(user.model),'Model must be a scalar struct.');
physicalDefaults=struct( ...
 'soilViscosity',struct('enabled',true,'referenceTemperature_C',20), ...
 'leakage',struct('jointStorage',struct('enabled',true,'tortuosity',1)));
model=joint_supply_th_model(combine(physicalDefaults,user.model));
if ~isfield(model,'dtMax'), model.dtMax=4*86400; end
assert(isscalar(model.dtMax)&&isfinite(model.dtMax)&& ...
 model.dtMax>0&&model.dtMax<=4*86400, ...
 'The regime-screening maximum step is four days.');
user.model=model;
if ~isfield(user,'outputFolder')
 user.outputFolder=fullfile(root,'segment_year20_joint_supply_8segment_illustrative');
end
if ~isfield(user,'maxYears'), user.maxYears=20; end
if ~isfield(user,'consecutiveYears'), user.consecutiveYears=21; end
if ~isfield(user,'productionYears'), user.productionYears=1; end
assert(user.maxYears<=20, ...
 'This age-20 wrapper permits at most 20 years; use the general runner otherwise.');
started=tic;
result=run_TH_segment_spinup_v14(user);
result.elapsedWallTime_s=toc(started);
result.jointOnlyModel=true;
result.jointFreezeModel=model.leakage.jointFreeze.enabled;
result.waterSupplyDiagnostics=true;
result.jointParameterStatus=['Eight-segment reference geometry; assumed defective ' ...
 'joint; reversible microchannel icing, finite channel storage and soil viscosity. ' ...
 'Soil, climate, groundwater and joint inputs are uncalibrated scenarios.'];
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
 assert(size(leak.hydraulics.history,1)==size(result.lastYear.history,1), ...
  'Annual water-supply diagnostics do not match the accepted-step history.');
 result.year20AnnualLeakage_L_per_m=1000*total;
 result.year20JointLeakage_L_per_m=1000*joint;
 result.year20BodyLeakage_L_per_m=1000*body;
 result.year20ByJoint_L_per_m=1000*leak.volumeByJoint_m3_per_m;
 result.year20ByMonth_L_per_m=1000*leak.monthlyTotalVolume_m3_per_m;
 result.year20MonthlyByJoint_L_per_m=1000*leak.monthlyVolumeByJoint_m3_per_m;
 result.year20DailyLeakage_L_per_day_per_m=1000*leak.dailyTotalVolume_m3_per_m;
end
% Keep wrapper metadata after the generic annual runner's completed-year save.
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
