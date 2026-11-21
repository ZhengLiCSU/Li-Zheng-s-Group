function zFull = lift_premelting_symmetry_state(z,model)
%LIFT_PREMELTING_SYMMETRY_STATE Expand a cyclic state to all eight segments.
arguments
    z (:,1) double {mustBeFinite}
    model (1,1) struct
end
assert(isfield(model,'symmetry'),'PhaseSwitch:Symmetry','A cyclic model is required.');
n=numel(model.weight); ns=numel(model.capacity);
assert(numel(z)==6*n+ns,'PhaseSwitch:Shape','Wrong cyclic state length.');
Q=model.symmetry.cellLift; Pt=model.symmetry.temperatureLift;
zFull=zeros(6*size(Q,1)+size(Pt,1),1);
nf=size(Q,1);
for k=0:5, zFull(k*nf+(1:nf))=Q*z(k*n+(1:n)); end
zFull(6*nf+(1:size(Pt,1)))=Pt*z(6*n+(1:ns));
end
