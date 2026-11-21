function s=thm_constitutive_signature(p)
%THM_CONSTITUTIVE_SIGNATURE Physical settings and nonlegacy closure guards.
% Defaults at a 0-C joint bulk freezing point retain the legacy signature:
% their equations are unchanged by the corrected option forwarding. A
% nonzero point requires the corrected common-point closure, preventing
% silent continuation of an older checkpoint that used 0 C for its slit.
nonphysical={'mesh','meshFile','durationDays','dtInitial','dtMax','dtMin', ...
    'outputEveryDays','outputPrefix','saveOutput','plot','verbose','initial', ...
    'solver','spinup'};
assert(isstruct(p)&&isscalar(p)&&all(isfield(p,nonphysical)), ...
    'THM:SignatureParameters','The complete resolved parameter structure is required.');
s=rmfield(p,nonphysical);
if isfield(s.leakage.jointFreeze,'flowWeightedRadiusProfile_m')
    if isempty(s.leakage.jointFreeze.flowWeightedRadiusProfile_m)
        s.leakage.jointFreeze=rmfield(s.leakage.jointFreeze,'flowWeightedRadiusProfile_m');
    else
        s.jointSpatialGeometryRevision='fixed_connected_area_radius_profile_20261007';
    end
end
% An absent spatial field preserves historical uniform-material signatures.
if isfield(s,'soilKsatField')&&isfield(s.soilKsatField,'quadratureCoordinates_m') ...
        &&isempty(s.soilKsatField.quadratureCoordinates_m)&&isempty(s.soilKsatField.quadratureLogMultiplier)
    s.soilKsatField=rmfield(s.soilKsatField,{'quadratureCoordinates_m','quadratureLogMultiplier'});
end
if isfield(s,'soilKsatField')&&isempty(s.soilKsatField.coordinates_m) ...
        &&isempty(s.soilKsatField.logMultiplier)
    s=rmfield(s,'soilKsatField');
end
% Disabled fixed-damage geometry leaves all historic equations unchanged.
if isfield(s.leakage,'jointDamage')
    if ~s.leakage.jointDamage.enabled
        s.leakage=rmfield(s.leakage,'jointDamage');
    else
        s.jointDamageGeometryClosureRevision='fixed_capillary_geometry_20261002';
    end
end
if isfield(s.leakage,'jointCycleDamage')
    if ~s.leakage.jointCycleDamage.enabled
        s.leakage=rmfield(s.leakage,'jointCycleDamage');
    else
        s.jointCycleDamageClosureRevision='constant_volume_local_ice_cycles_20261002';
    end
end
if p.thm.enabled&&p.leakage.jointFreeze.bulkFreezingTemperature_C~=0
    s.jointSlitFreezingClosureRevision='common_bulk_freezing_point_20261001';
end
end
