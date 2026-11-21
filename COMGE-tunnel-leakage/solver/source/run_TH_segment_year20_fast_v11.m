function result=run_TH_segment_year20_fast_v11(user)
%RUN_TH_SEGMENT_YEAR20_FAST_V11 Age-specific 20-year annual-leakage sample.
% The 4-day maximum is a screening configuration validated against the
% supplied deterministic v10 example; revalidate it for other parameters.
% Supply a unique outputFolder for each sample and use resume=true after
% interruption. The ordinary v11 annual runner retains a 1-day default.
if nargin<1, user=struct(); end
assert(isstruct(user) && isscalar(user));
root=fileparts(mfilename('fullpath'));
if ~isfield(user,'outputFolder')
    user.outputFolder=fullfile(root,'segment_year20_fast_v11_output');
end
if ~isfield(user,'model'), user.model=struct(); end
assert(isstruct(user.model) && isscalar(user.model));
if ~isfield(user.model,'dtMax'), user.model.dtMax=4*86400; end
assert(user.model.dtMax>0 && user.model.dtMax<=4*86400);
user.maxYears=20;
user.consecutiveYears=21; % Do not mistake age-20 output for periodic state.
user.productionYears=1;
result=run_TH_segment_spinup_v11(user);
if result.completedYears==20
    result.year20AnnualLeakage_L_per_m= ...
        1000*result.lastYear.leakage.totalVolume_m3_per_m;
else
    result.year20AnnualLeakage_L_per_m=NaN;
end
end
