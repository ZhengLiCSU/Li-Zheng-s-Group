function report = audit_joint_channel_geometry(inventory, jointWidth, liningThickness)
%AUDIT_JOINT_CHANNEL_GEOMETRY Bound equivalent water inventory by its band.
% Each longitudinal drain is represented by a finite circumferential band.
% Its connected void volume per metre of tunnel axis cannot exceed width*L.
% This checks the model's equivalent geometry, not measured gasket porosity.
arguments
    inventory (1,1) struct
    jointWidth double {mustBeReal, mustBeFinite, mustBePositive}
    liningThickness double {mustBeReal, mustBeFinite, mustBePositive}
end
assert(isfield(inventory,'waterVolume_m2'), 'JointGeometry:Inventory', ...
    'Inventory must contain connected water volume per metre of tunnel axis.');
volume = inventory.waterVolume_m2(:);
validateattributes(volume, {'numeric'}, {'real','finite','nonnegative'});
numberOfJoints = numel(volume);
width = jointVector(jointWidth, numberOfJoints);
thickness = jointVector(liningThickness, numberOfJoints);
availableVolume = width.*thickness;
fraction = volume./availableVolume;
assert(all(fraction <= 1 + 1e-10), 'JointGeometry:ExcessVolume', ...
    ['Connected channel volume exceeds its finite joint band. ' ...
     'Reconsider joint conductance, radius distribution and tortuosity together.']);
report = struct('passed',true,'waterVolume_m2',volume, ...
    'availableBandVolume_m2',availableVolume,'equivalentVoidFraction',fraction, ...
    'scope','Model-equivalent channel inventory; finite band footprint; no gasket porosity measurement.');
end

function value = jointVector(value, count)
value = value(:);
if isscalar(value)
    value = repmat(value, count, 1);
end
assert(numel(value) == count, 'JointGeometry:Shape', ...
    'Supply one width/thickness or one value per joint.');
end
