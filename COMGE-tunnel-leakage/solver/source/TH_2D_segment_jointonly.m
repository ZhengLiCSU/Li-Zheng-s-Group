function out=TH_2D_segment_jointonly(opt)
%TH_2D_SEGMENT_JOINTONLY Two-way THM with impermeable concrete and eight joints.
% Uses v16: conservative SCV/CVFE TH, poroelastic soil, eight independently
% meshed contact segments, bolts and variable-aperture joint hydraulic paths.
% Default joint values are illustrative and are not a tunnel prediction.
if nargin<1, opt=struct(); end
if ischar(opt) && strcmp(opt,'selftest')
    out=TH_2D_segment_v16('thmselftest');
    return
end
opt=joint_contact_thm_model(opt);
out=TH_2D_segment_v16(opt);
assert(out.parameters.liningK==0 && ...
    out.leakage.bodyVolume_m3_per_m==0, ...
    'Concrete body leakage was detected in the joint-only model.');
out.jointOnlyModel=true;
end
