function result=run_TH_segment_year20_cvfe(user)
%RUN_TH_SEGMENT_YEAR20_CVFE Full age history with phase-integrated CVFE.
% Separate schema/output from archived v14. Default maximum step one day.
% This remains equal-density TH, with no deformation or cyclic damage law.
if nargin<1, user=struct(); end
assert(isstruct(user)&&isscalar(user));
root=fileparts(mfilename('fullpath'));
if ~isfield(user,'model'), user.model=struct(); end
model=joint_supply_th_model(user.model);
if ~isfield(model,'dtMax'), model.dtMax=86400; end
if ~isfield(model,'integration'), model.integration=struct('soilQuadratureOrder',5); end
user.model=model;
if ~isfield(user,'outputFolder')
    user.outputFolder=fullfile(root,'segment_year20_v15_cvfe_illustrative');
end
if ~isfield(user,'maxYears'), user.maxYears=20; end
if ~isfield(user,'consecutiveYears'), user.consecutiveYears=21; end
if ~isfield(user,'productionYears'), user.productionYears=1; end
assert(user.maxYears<=20);
started=tic; result=run_TH_segment_spinup_v15(user);
result.elapsedWallTime_s=toc(started);
result.jointOnlyModel=true;
result.modelStatus='Conservative equal-density TH; example inputs; mechanical feedback and gasket calibration absent.';
result.year20AnnualLeakage_L_per_m=NaN;
if result.completedYears==20
    result.year20AnnualLeakage_L_per_m=1000*result.lastYear.leakage.totalVolume_m3_per_m;
end
assert(result.lastYear.parameters.liningK==0&&result.lastYear.leakage.bodyVolume_m3_per_m==0);
save(fullfile(user.outputFolder,'TH_segment_spinup_results.mat'),'result','-v7');
end
