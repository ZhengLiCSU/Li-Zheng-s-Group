function [v,schema]=frozen_soil_jit_parameters(p)
%FROZEN_SOIL_JIT_PARAMETERS Pack scientific scalars for the isolated JIT law.
% No mesh, timestep or solver setting is encoded. This packer performs the
% existing public law's validation once, before repeated bulk evaluations.
% Reference enthalpy and disabled/missing viscosity use the original defaults.

names={'rho','g','L','Tm','porosity','thetaR','vgAlpha','vgN','Ss', ...
    'Ksat','iceImpedance','cw','ci','rhoS','cs','kw','ki','ks','ka'};
assert(isstruct(p)&&isscalar(p),'FrozenSoilJIT:Parameters', ...
    'Supply a scalar physical-parameter structure.');
% Reuse its public validation and avoid introducing a weaker material domain.
frozen_soil_constitutive(4,0,p,false);
v=zeros(22,1);
for j=1:numel(names), v(j)=p.(names{j}); end
v(20)=0; v(21)=20; v(22)=0;
if isfield(p,'soilViscosity')&&p.soilViscosity.enabled
    v(20)=1; v(21)=p.soilViscosity.referenceTemperature_C;
end
if isfield(p,'enthalpyReference_Jkg'), v(22)=p.enthalpyReference_Jkg; end
schema.version='frozen_soil_properties_v1_20260930';
schema.parameterNames=[names,{'viscosityEnabled', ...
    'viscosityReferenceTemperature_C','enthalpyReference_Jkg'}];
schema.propertyNames={'sl','si','sg','thetaW','thetaI','thetaA', ...
    'total','elasticWater','mass','energy','k','K'};
schema.parameterCount=22; schema.propertyCount=12;
schema.note=['Properties-only equal-density frozen law. Matrix rows are ', ...
    'constitutive points; SCV interpolation and integration remain external.'];
end
