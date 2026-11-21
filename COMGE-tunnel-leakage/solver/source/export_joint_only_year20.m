function files=export_joint_only_year20(outputFolder)
%EXPORT_JOINT_ONLY_YEAR20 Write reviewable age-20 joint leakage artifacts.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(root,'segment_year20_jointonly_8segment_illustrative');
end
s=load(fullfile(outputFolder,'TH_segment_spinup_results.mat'),'result');
r=s.result;
assert(r.completedYears==20 && r.lastYear.parameters.liningK==0);
leak=r.lastYear.leakage;
assert(leak.bodyVolume_m3_per_m==0 && ...
    abs(sum(leak.volumeByJoint_m3_per_m)-leak.totalVolume_m3_per_m)<1e-9);
jointNames=arrayfun(@(k)sprintf('Joint%d_L_per_m',k), ...
    1:numel(leak.volumeByJoint_m3_per_m),'UniformOutput',false);
monthly=array2table([(1:12)' 1000*leak.monthlyTotalVolume_m3_per_m(:) ...
    1000*leak.monthlyVolumeByJoint_m3_per_m], ...
    'VariableNames',[{'Month','Total_L_per_m'},jointNames]);
nj=numel(leak.volumeByJoint_m3_per_m);
config=r.lastYear.parameters.leakage;
angles=config.jointAnglesRad(:);
conductance=config.jointConductance(:);
if isscalar(conductance), conductance=repmat(conductance,nj,1); end
entryHead=config.jointEntryHead(:);
if isempty(entryHead), entryHead=zeros(nj,1);
elseif isscalar(entryHead), entryHead=repmat(entryHead,nj,1); end
annual=table((1:nj)',angles,angles*180/pi,conductance,entryHead, ...
    1000*leak.volumeByJoint_m3_per_m(:), ...
    'VariableNames',{'Joint','Angle_rad','Angle_deg', ...
    'Conductance_m_per_s','EntryHead_m','Volume_L_per_m'});
files=struct('monthly',fullfile(outputFolder,'year20_monthly_joint_leakage.csv'), ...
    'annual',fullfile(outputFolder,'year20_annual_joint_leakage.csv'), ...
    'figure',fullfile(outputFolder,'year20_joint_leakage.png'));
writetable(monthly,files.monthly);
writetable(annual,files.annual);
fig=plot_joint_only_year20(r,files.figure);
close(fig);
end
