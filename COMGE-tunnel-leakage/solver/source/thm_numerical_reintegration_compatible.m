function compatible=thm_numerical_reintegration_compatible(old,new)
%THM_NUMERICAL_REINTEGRATION_COMPATIBLE Only phase/path quadrature may change.
% Both inputs are physical signatures. Soil/boundary/contact/closure fields
% remain exact. This is for a NEW conditional run from a saved complete X,
% never for relabelling or recovering an existing accepted-step inventory.
old=removeNumerical(old);new=removeNumerical(new);
compatible=isequaln(old,new);
end
function s=removeNumerical(s)
assert(isstruct(s)&&isscalar(s)&&isfield(s,'integration') ...
    &&isfield(s,'leakage')&&isfield(s.leakage,'jointFreeze'), ...
    'THM:ReintegrationSignature','Complete signatures are required.');
s=rmfield(s,'integration');
s.leakage.jointFreeze=rmfield(s.leakage.jointFreeze,'quadraturePoints');
end
