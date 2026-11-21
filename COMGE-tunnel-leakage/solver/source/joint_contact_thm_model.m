function opt=joint_contact_thm_model(user)
%JOINT_CONTACT_THM_MODEL Example eight-joint monolithic THM configuration.
% Retains the refined joint-only supply geometry and its stated parameters.
% Mechanical/contact/aperture defaults reside in TH_2D_segment_v16. Override
% them through user.thm / user.mechanics after calibration. Historical v14
% TH remains available explicitly as TH_2D_segment_v14(joint_supply_th_model).
if nargin<1, user=struct(); end
opt=joint_supply_th_model(user);
if ~isfield(opt,'thm'), opt.thm=struct(); end
if ~isfield(opt.thm,'enabled'), opt.thm.enabled=true; end
assert(opt.thm.enabled,'JointContactTHM:Model', ...
    'The current entry requires two-way THM. Use the explicit v15/v14 core for TH comparisons.');
opt.liningK=0;
end
