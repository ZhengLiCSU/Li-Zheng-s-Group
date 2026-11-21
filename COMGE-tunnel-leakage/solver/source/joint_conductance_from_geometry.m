function [conductance, report] = joint_conductance_from_geometry(geometry, freeze)
%JOINT_CONDUCTANCE_FROM_GEOMETRY Generate admissible correlated joint inputs.
% Geometry: bandWidth_m, liningThickness_m, equivalentVoidFraction [0,1],
% tortuosity>=1, density_kg_m3 and gravity_m_s2 (positive scalars).
% Freeze supplies the r^4-weighted lognormal radius, width and mu reference.
% Specified void fraction fixes V=phi*width*L; C scales as rQ^2/tau^2.
% Thus random radius/void fraction inputs need not invent an independent
% conductance and inconsistent water inventory. These are equivalent network
% coordinates, not molecular porosity of solid EPDM. No fitting is done.
arguments
    geometry (1,1) struct
    freeze (1,1) struct
end
keys = {'bandWidth_m','liningThickness_m','equivalentVoidFraction', ...
    'tortuosity','density_kg_m3','gravity_m_s2'};
assert(all(isfield(geometry,keys)), 'JointGeometry:Inputs','Supply all six unit-labelled geometry fields.');
freeze = joint_freeze_defaults(freeze);
arrays = {geometry.bandWidth_m,geometry.liningThickness_m, ...
    geometry.equivalentVoidFraction,geometry.tortuosity, ...
    freeze.flowWeightedRadius_m,freeze.logRadiusStd};
count = max(cellfun(@numel,arrays));
resolved = cell(size(arrays));
for index = 1:numel(arrays)
    value = arrays{index};
    validateattributes(value,{'numeric'},{'vector','real','finite','nonempty'});
    value = value(:);
    if isscalar(value)
        value = repmat(value,count,1);
    end
    assert(numel(value)==count,'JointGeometry:Shape','Geometry distributions must have matching dimensions.');
    resolved{index} = value;
end
width = resolved{1};
thickness = resolved{2};
fraction = resolved{3};
tortuosity = resolved{4};
radius = resolved{5};
sigma = resolved{6};
assert(all(fraction>=0 & fraction<=1),'JointGeometry:VoidFraction','Equivalent void fraction must lie in [0,1].');
assert(all(width>0 & thickness>0 & radius>0 & sigma>0 & tortuosity>=1), ...
    'JointGeometry:Range','Positive dimensions/radius/width and tortuosity>=1 are required.');
validateattributes(geometry.density_kg_m3,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(geometry.gravity_m_s2,{'numeric'},{'real','finite','scalar','positive'});
mu = water_viscosity_iapws(freeze.referenceTemperature_C);
conductance = geometry.density_kg_m3*geometry.gravity_m_s2 ...
    *fraction.*width.*radius.^2.*exp(-2*sigma.^2)./(8*mu*thickness.*tortuosity.^2);
assert(all(isfinite(conductance)), 'JointGeometry:Overflow','Geometry produces nonfinite conductance.');
volume = fraction.*width.*thickness;
report = audit_joint_channel_geometry(struct('waterVolume_m2',volume),width,thickness);
report.referenceConductance_m_s = conductance;
report.connectedArea_m2_per_m = volume./(thickness.*tortuosity);
report.scope = 'Equivalent lognormal channel network generated from bounded void inventory; correlated C/radius/tortuosity; no empirical gasket calibration.';
end
