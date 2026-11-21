function [pField,metadata] = apply_premelting_aperture_field(model,p,field)
%APPLY_PREMELTING_APERTURE_FIELD Mesh-independent, positive geometric input.
% field.xi ranges from outer face 0 to inner face 1. aperture_m is 8-by-K.
% Log-linear interpolation samples the SAME material realization on any
% mesh. This supplies geometry; it does not assign a probability model.
arguments
    model (1,1) struct
    p (1,1) struct
    field (1,1) struct
end
assert(~isfield(model,'symmetry'),'PhaseSwitch:Field','Map the field on the full eight-joint mesh before reducing a nominal case.');
assert(isequal(sort(fieldnames(field)),sort({'xi';'aperture_m'})),'PhaseSwitch:Field','Supply only xi and aperture_m.');
xi=field.xi(:)'; values=field.aperture_m;
assert(isreal(xi) && numel(xi)>=2 && all(isfinite(xi)) && xi(1)==0 && xi(end)==1 && all(diff(xi)>0), ...
    'PhaseSwitch:Field','Strictly increasing material coordinates must span [0,1].');
assert(isreal(values) && isequal(size(values),[8,numel(xi)]) && all(isfinite(values)&values>0,'all'), ...
    'PhaseSwitch:Field','Supply positive finite aperture values on eight material paths.');
nq=model.nq; query=((1:nq)-.5)/nq;
sample=zeros(nq,8);
for j=1:8, sample(:,j)=exp(interp1(xi,log(values(j,:)),query,'linear'))'; end
% Preserve constant inputs exactly for zero-variance regression.
for j=1:8
    if all(values(j,:)==values(j,1)), sample(:,j)=values(j,1); end
end
pField=p; pField.aperture0_m=sample(:);
metadata=struct('definition',field,'cellXi',repmat(query',8,1), ...
    'cellJointID',repelem((1:8)',nq),'cellAperture_m',sample(:), ...
    'scope','Deterministic realization of residual hydraulic aperture, sampled in material coordinates; no independent conductance or storage multiplier.');
end
