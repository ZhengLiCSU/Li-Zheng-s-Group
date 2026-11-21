function model=joint_only_th_model(overrides)
%JOINT_ONLY_TH_MODEL Example 2D TH case with impermeable concrete segments.
% Water enters the tunnel only through longitudinal segment joints. Each
% joint is an equivalent, finite-support interface drain per metre of
% tunnel axis, not a resolved gasket or a circumferential (ring) joint.
%
% The default 5+2+1 angles follow the Shiziyang prototype drawing, scaled
% onto this model's different 12-m-OD illustrative mesh. Only ONE joint is
% assigned an assumed microchannel defect; the other seven are sealed.
% This is NOT a calibrated prediction for Shiziyang or a cold-region site.
% Joints extend only one ring width; staggered rings remain unresolved.
% Hydraulic parameters below must be calibrated for a real tunnel.
% jointWidth_m locates
% the drainage band on the circumference; it is NOT joint opening.
% jointConductance [m/s] is the integrated conductance of one joint per
% metre of tunnel axis: q_j [m^2/s] = C_j [m/s] * max(dH-Hentry,0) [m].
if nargin<1, overrides=struct(); end
assert(isstruct(overrides)&&isscalar(overrides),'Overrides must be a scalar struct.');
layout=hsr_eight_segment_layout();
defectConductance=zeros(1,8);
defectConductance(strcmp(layout.jointLabels,'B5-B4'))=2e-10;
% The silt retention shape follows Dall'Amico et al. (2011), Table 1;
% porosity 0.35 is a hypothetical denser/compacted soil scenario, not that
% paper's 0.49 material value. Both require site-specific measurements.
% https://tc.copernicus.org/articles/5/469/2011/tc-5-469-2011.pdf
% Ksat, air temperatures and 5-m water table are explicitly hypothetical
% cold-region scenario choices, not site measurements. The 0.06 W/m^2
% basal flux represents the Chinese continental non-anomalous mean;
% the initial/deep-inflow gradient uses the same unfrozen bulk k.
soilBulkK=(1-.35)*2.3+.35*.57;
referenceGradient=.06/soilBulkK;
model=struct('liningK',0,'porosity',.35,'thetaR',.05, ...
    'vgAlpha',.65,'vgN',1.67,'Ksat',1e-7, ...
    'geothermalFlux',.06,'waterTableDepth',5, ...
    'topMin',-10,'topMax',15,'tunnelMin',-5,'tunnelMax',8, ...
    'groundwaterInflowT',2.5,'groundwaterInflowGradient',referenceGradient, ...
    'initial',struct('surfaceT',2.5,'gradient',referenceGradient), ...
    'leakage',struct('segmentCount',8,'segmentStartAngleRad',0, ...
    'interiorHead',0,'jointAnglesRad',layout.jointAnglesRad, ...
    'jointWidth_m',.02,'jointConductance',defectConductance, ...
    'jointEntryHead',0));
model=combine(model,overrides);
if isfield(overrides,'leakage') && ...
        (isfield(overrides.leakage,'segmentCount') || ...
         isfield(overrides.leakage,'segmentStartAngleRad')) && ...
        ~isfield(overrides.leakage,'jointAnglesRad')
    n=model.leakage.segmentCount;
    if n==8
        layout=hsr_eight_segment_layout(model.leakage.segmentStartAngleRad);
        model.leakage.jointAnglesRad=layout.jointAnglesRad;
        if ~isfield(overrides.leakage,'jointConductance')
            defectConductance=zeros(1,8);
            defectConductance(strcmp(layout.jointLabels,'B5-B4'))=2e-10;
            model.leakage.jointConductance=defectConductance;
        end
    else
        model.leakage.jointAnglesRad=model.leakage.segmentStartAngleRad+ ...
            (0:n-1)*2*pi/n;
        if ~isfield(overrides.leakage,'jointConductance')
            model.leakage.jointConductance=0;
        end
    end
end
assert(isscalar(model.liningK)&&isfinite(model.liningK)&&model.liningK==0, ...
    'Joint-only model requires liningK=0: concrete body is impermeable.');
L=model.leakage;
assert(~isempty(L.jointAnglesRad),'Specify at least one longitudinal joint.');
assert(all(isfinite(L.jointAnglesRad(:))),'Joint angles must be finite.');
assert(all(isfinite(L.jointWidth_m(:)))&&all(L.jointWidth_m(:)>0), ...
    'Joint support width must be positive.');
assert(all(isfinite(L.jointConductance(:)))&&all(L.jointConductance(:)>=0), ...
    'Joint conductance must be nonnegative.');
assert(all(isfinite(L.jointEntryHead(:)))&&all(L.jointEntryHead(:)>=0), ...
    'Joint entry head must be nonnegative.');
end

function a=combine(a,b)
fields=fieldnames(b);
for k=1:numel(fields)
    key=fields{k};
    if isfield(a,key)&&isstruct(a.(key))&&isstruct(b.(key))
        a.(key)=combine(a.(key),b.(key));
    else
        a.(key)=b.(key);
    end
end
end
