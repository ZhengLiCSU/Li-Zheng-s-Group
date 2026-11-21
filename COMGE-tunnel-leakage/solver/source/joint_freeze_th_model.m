function model=joint_freeze_th_model(overrides)
%JOINT_FREEZE_TH_MODEL Current joint-only model with reversible freeze flow.
% The hydraulic conductance and microchannel radius distribution remain
% illustrative inputs requiring independent joint tests for calibration.
if nargin<1, overrides=struct(); end
model=joint_only_th_model(overrides);
settings=joint_freeze_defaults(struct('enabled',true));
if isfield(model.leakage,'jointFreeze')
    names=fieldnames(model.leakage.jointFreeze);
    for k=1:numel(names)
        assert(isfield(settings,names{k}),'Unknown joint-freezing setting.');
        settings.(names{k})=model.leakage.jointFreeze.(names{k});
    end
end
model.leakage.jointFreeze=joint_freeze_defaults(settings);
end
