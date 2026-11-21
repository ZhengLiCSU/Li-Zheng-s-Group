function result=run_TH_segment_year20_jointonly_v11(user)
%RUN_TH_SEGMENT_YEAR20_JOINTONLY_V11 Legacy constant-conductance joint baseline.
% The concrete body is hydraulically impermeable. The default is an
% eight-segment, one-defective-joint hypothetical cold-silt scenario.
% Replace all site-dependent soil, climate, water and joint fields before
% interpreting leakage as an engineering prediction.
% The inherited 4-day maximum step is a screening setting; verify it for
% the chosen joint case against a smaller time step.
% Example: r=run_TH_segment_year20_jointonly_v11;
% Resume:  r=run_TH_segment_year20_jointonly_v11(struct('resume',true));
if nargin<1, user=struct(); end
assert(isstruct(user)&&isscalar(user),'Input must be a scalar struct.');
root=fileparts(mfilename('fullpath'));
if ~isfield(user,'model'), user.model=struct(); end
user.model=joint_only_th_model(user.model);
if ~isfield(user,'outputFolder')
    user.outputFolder=fullfile(root,'segment_year20_jointonly_8segment_illustrative');
end
result=run_TH_segment_year20_fast_v11(user);
result.jointOnlyModel=true;
result.jointParameterStatus=['Eight-segment reference geometry; one assumed ' ...
    'defective joint; soil, climate and water boundaries are uncalibrated.'];
result.year20JointLeakage_L_per_m=NaN;
result.year20BodyLeakage_L_per_m=NaN;
result.year20ByJoint_L_per_m=NaN;
result.year20ByMonth_L_per_m=NaN;
if result.completedYears==20
    leak=result.lastYear.leakage;
    total=leak.totalVolume_m3_per_m;
    joint=sum(leak.volumeByJoint_m3_per_m);
    body=leak.bodyVolume_m3_per_m;
    assert(result.lastYear.parameters.liningK==0 && ...
        abs(body)<1e-12 && abs(joint-total)<1e-9, ...
        'Joint-only mass partition failed.');
    assert(abs(sum(leak.volumeBySegment_m3_per_m)-total)<1e-9, ...
        'Segment rates do not partition total joint leakage.');
    assert(abs(sum(leak.monthlyTotalVolume_m3_per_m)-total)<1e-9, ...
        'Monthly rates do not partition total joint leakage.');
    result.year20JointLeakage_L_per_m=1000*joint;
    result.year20BodyLeakage_L_per_m=1000*body;
    result.year20ByJoint_L_per_m=1000*leak.volumeByJoint_m3_per_m;
    result.year20ByMonth_L_per_m=1000*leak.monthlyTotalVolume_m3_per_m;
end
% The generic annual runner writes this file before the joint-only fields
% are added. Save the final result so reloading preserves the same metadata.
save(fullfile(user.outputFolder,'TH_segment_spinup_results.mat'),'result','-v7');
end
