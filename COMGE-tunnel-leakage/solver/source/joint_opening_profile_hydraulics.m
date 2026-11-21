function [C,dCdg,dCdT,D]=joint_opening_profile_hydraulics(Cref,gapProfile,Tprofile,options)
%JOINT_OPENING_PROFILE_HYDRAULICS Variable-aperture slit series resistance.
% [C,dCdg,dCdT,D]=joint_opening_profile_hydraulics(Cref,gapProfile,Tprofile,options)
% accepts nj-by-nq signed mechanical openings and temperatures along each
% joint thickness. Physical options and local slit thermodynamics are those
% of joint_opening_hydraulics; the original capillary defect remains a
% separate parallel branch Cref*baseMobilityFactor. C is nj-by-1 [m/s];
% dCdg and dCdT are nj-by-nq partial derivatives at the supplied samples.
%
% Every local slit conductance c_k uses the FULL path length L, then
% Copen=1/sum_k(w_k/c_k), where positive weights are length fractions.
% This is the quadrature of local cubic-law hydraulic resistance; no minimum
% aperture, averaged aperture, artificial conductance floor or pre-clipping
% of a positive opening is substituted. Any genuinely closed/frozen local
% passage makes the open branch exactly zero; the base defect is unchanged.
% dCopen/dc_k=w_k*(Copen/c_k)^2 on an unblocked branch.
%
% Cref is scalar or nj values; gapProfile and Tprofile are the same shape.
% options.baseMobilityFactor is scalar or nj values; its optional temperature
% derivative is nj-by-nq. options.sampleWeights is a positive nq vector.
% Per-joint pathLength_m and openingTransmissionFraction remain supported.
% options.openPathMask (scalar or nj logical values, default true) supplies
% an independent whole-path continuity check, e.g. all mechanical radial
% nodal gaps>0. A false mask blocks only new hydraulic flow and its tangents;
% local apertures and partially open channel inventories remain unchanged.
% Temperature may be [] only when freezing and viscosity are both disabled.
%
% V=L*sum_k(w_k*b_k), sample liquid/ice volumes and their derivatives are
% returned for a separate channel-inventory equation. This module alone
% does not close moving-channel pressure/storage, inlet reversal or heat
% conservation. The slit thermodynamic regularization and aperture mapping
% retain the calibration limitations documented by the scalar helper.

if nargin<3, Tprofile=[]; end
if nargin<4 || isempty(options), options=struct(); end
assert(isstruct(options)&&isscalar(options),'JointOpeningProfile:Options', ...
    'options must be a scalar struct.');
assert(isnumeric(gapProfile)&&isreal(gapProfile)&&ismatrix(gapProfile) ...
    &&~isempty(gapProfile)&&all(isfinite(gapProfile(:))), ...
    'JointOpeningProfile:Gap','gapProfile must be a finite nonempty joint-by-sample matrix.');
[nj,nq]=size(gapProfile); gapProfile=double(gapProfile);
Cref=jointValue(Cref,nj,'Cref',0,Inf);
thermalProvided=~isempty(Tprofile);
if ~thermalProvided
    freeze=isfield(options,'freezeEnabled')&&options.freezeEnabled;
    viscosity=isfield(options,'viscosityEnabled')&&options.viscosityEnabled;
    assert(~freeze&&~viscosity,'JointOpeningProfile:MissingTemperature', ...
        'Enabled freezing or viscosity requires a temperature profile.');
    Tprofile=zeros(nj,nq);
end
assert(isnumeric(Tprofile)&&isreal(Tprofile)&&isequal(size(Tprofile),[nj nq]) ...
    &&all(isfinite(Tprofile(:))),'JointOpeningProfile:Temperature', ...
    'Tprofile must match the finite gapProfile matrix.');
Tprofile=double(Tprofile);
weights=ones(1,nq)/nq;
if isfield(options,'sampleWeights')&&~isempty(options.sampleWeights)
    value=options.sampleWeights;
    assert(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==nq ...
        &&all(isfinite(value(:)))&&all(value(:)>0),'JointOpeningProfile:Weights', ...
        'sampleWeights must contain positive finite length fractions.');
    weights=double(value(:)'); weights=weights/sum(weights);
end
baseFactor=ones(nj,1);
if isfield(options,'baseMobilityFactor')
    baseFactor=jointValue(options.baseMobilityFactor,nj,'baseMobilityFactor',0,Inf);
end
baseDerivative=zeros(nj,nq);
if isfield(options,'baseMobilityDerivative_dT')&&~isempty(options.baseMobilityDerivative_dT)
    baseDerivative=options.baseMobilityDerivative_dT;
    assert(isnumeric(baseDerivative)&&isreal(baseDerivative) ...
        &&isequal(size(baseDerivative),[nj nq])&&all(isfinite(baseDerivative(:))), ...
        'JointOpeningProfile:BaseDerivative','Base derivative must be nj-by-nq.');
    baseDerivative=double(baseDerivative);
end
% The scalar helper evaluates a single local aperture/temperature for every
% flattened point. Its local path length stays L, not w_k*L: weights below
% assemble resistance lengths once. MATLAB column-major flattening is used.
localOptions=options;
openPathMask=true(nj,1);
if isfield(localOptions,'openPathMask')
    value=localOptions.openPathMask;
    assert((islogical(value)||isnumeric(value))&&isreal(value)&&isvector(value) ...
        &&all(isfinite(value(:)))&&all(ismember(value(:),[0 1])) ...
        &&(isscalar(value)||numel(value)==nj),'JointOpeningProfile:OpenPathMask', ...
        'openPathMask must be logical scalar or one logical value per joint.');
    openPathMask=logical(value(:));
    if isscalar(openPathMask), openPathMask=repmat(openPathMask,nj,1); end
    localOptions=rmfield(localOptions,'openPathMask');
end
inventoryNames={'pathLength_m','openingTransmissionFraction'};
for k=1:numel(inventoryNames)
    name=inventoryNames{k};
    if isfield(localOptions,name)
        minimum=0; maximum=1;
        if strcmp(name,'pathLength_m'), minimum=realmin; maximum=Inf; end
        value=jointValue(localOptions.(name),nj,name,minimum,maximum);
        localOptions.(name)=repmat(value,nq,1);
    end
end
localOptions.sampleWeights=1;
localOptions.baseMobilityFactor=0;
localOptions.baseMobilityDerivative_dT=zeros(nj*nq,1);
[localC,localGapDerivative,localTemperatureDerivative,localD]= ...
    joint_opening_hydraulics(repmat(Cref,nq,1),gapProfile(:),Tprofile(:),localOptions);
localC=reshape(localC,nj,nq);
localGapDerivative=reshape(localGapDerivative,nj,nq);
localTemperatureDerivative=reshape(localTemperatureDerivative,nj,nq);
[openC,chain]=seriesConductance(localC,weights);
openC(~openPathMask)=0; chain(~openPathMask,:)=0;
referenceLocal=reshape(localD.referenceOpenConductance,nj,nq);
referenceOpen=seriesConductance(referenceLocal,weights);
dCdg=chain.*localGapDerivative;
dCdT=chain.*localTemperatureDerivative+Cref.*baseDerivative;
C=Cref.*baseFactor+openC;
aperture=reshape(localD.aperture_m,nj,nq);
sampleVolume=reshape(localD.newslitGeometricVolume_m2,nj,nq).*weights;
dVolumeG=reshape(localD.dVolume_dGap_m,nj,nq).*weights;
phase=reshape(localD.liquidFractionBySample,nj,nq);
liquidVolume=reshape(localD.liquidVolumeBySample_m2,nj,nq).*weights;
iceVolume=reshape(localD.iceVolumeBySample_m2,nj,nq).*weights;
dLiquidVolumeG=reshape(localD.dLiquidVolume_dGap_m,nj,nq).*weights;
dIceVolumeG=reshape(localD.dIceVolume_dGap_m,nj,nq).*weights;
dLiquidVolumeT=reshape(localD.dLiquidVolume_dTemperature_m2_K,nj,nq).*weights;
factor=zeros(nj,1); openReference=referenceOpen>0;
factor(openReference)=openC(openReference)./referenceOpen(openReference);
baseAperture=reshape(localD.baseAperture,nj,nq);
status=repmat({'closed_or_disabled'},nj,1); status(Cref==0)={'sealed'};
status(openC>0)={'open_slit_profile'};
positivePath=all(aperture>0,2);
status(Cref>0&positivePath&openC==0)={'open_slit_ice_blocked'};
status(Cref>0&~openPathMask)={'whole_path_closed'};
s=localD.settings;
s.pathLength_m=reshape(jointValue(s.pathLength_m,nj*nq,'pathLength_m',realmin,Inf),nj,nq);
s.pathLength_m=s.pathLength_m(:,1);
s.openingTransmissionFraction=reshape(jointValue(s.openingTransmissionFraction, ...
    nj*nq,'openingTransmissionFraction',0,1),nj,nq);
s.openingTransmissionFraction=s.openingTransmissionFraction(:,1);
s.sampleWeights=weights; s.baseMobilityFactor=baseFactor;
s.baseMobilityDerivative_dT=baseDerivative;
s.openPathMask=openPathMask;
assumptions=localD.assumptions;
assumptions{2}='Unit axial width; laminar local cubic flow; variable aperture and path-length quadrature assembled as series resistance.';
assumptions{end+1}='Pointwise gap and temperature interpolation is supplied by the caller; this kernel does not collapse a profile to its minimum or mean.';
D=struct('conductance',C,'baseConductance',Cref.*baseFactor, ...
    'openConductance',openC,'referenceOpenConductance',referenceOpen, ...
    'openMobilityFactor',factor,'localConductance',localC, ...
    'localReferenceConductance',referenceLocal,'dSeries_dLocalConductance',chain, ...
    'baseAperture',baseAperture(:,1),'openingContribution',aperture,'aperture_m',aperture, ...
    'freezeTemperature_C',reshape(localD.freezeTemperature_C,nj,nq), ...
    'conductingFractionBySample',phase,'liquidFractionBySample',phase, ...
    'dLiquidFraction_dTemperature',reshape(localD.dLiquidFraction_dTemperature,nj,nq), ...
    'dLiquidFraction_dGap',reshape(localD.dLiquidFraction_dGap,nj,nq), ...
    'localMobility',reshape(localD.localMobility,nj,nq), ...
    'newslitGeometricVolume_m2',sum(sampleVolume,2), ...
    'newslitWaterVolume_m2',sum(sampleVolume,2),'dVolume_dGap_m',dVolumeG, ...
    'sampleVolume_m2',sampleVolume,'liquidVolumeBySample_m2',liquidVolume, ...
    'iceVolumeBySample_m2',iceVolume,'dLiquidVolume_dGap_m',dLiquidVolumeG, ...
    'dIceVolume_dGap_m',dIceVolumeG, ...
    'dLiquidVolume_dTemperature_m2_K',dLiquidVolumeT, ...
    'dIceVolume_dTemperature_m2_K',-dLiquidVolumeT, ...
    'sealedMask',Cref==0,'status',{status},'sampleWeights',weights, ...
    'openPathMask',openPathMask, ...
    'settings',s,'thermalProfileProvided',thermalProvided,'assumptions',{assumptions}, ...
    'baseApertureInterpretation',localD.baseApertureInterpretation, ...
    'hydraulicLawDOI',localD.hydraulicLawDOI, ...
    'slitThermodynamicsDOI',localD.slitThermodynamicsDOI, ...
    'constitutiveHelper','joint_opening_hydraulics');
end

function [C,derivative]=seriesConductance(localC,weights)
% Scale by the minimum positive value, retaining genuinely zero passages.
% C/localC is formed before squaring, avoiding underflow in C^2/localC^2.
nj=size(localC,1); nq=size(localC,2);
C=zeros(nj,1); derivative=zeros(nj,nq);
active=all(localC>0,2);
if any(active)
    values=localC(active,:); scale=min(values,[],2);
    scaledResistance=sum(weights.*(scale./values),2);
    C(active)=scale./scaledResistance;
    derivative(active,:)=weights.*(C(active)./values).^2;
end
end

function value=jointValue(value,nj,name,minimum,maximum)
assert(isnumeric(value)&&isreal(value)&&isvector(value)&&all(isfinite(value(:))) ...
    &&all(value(:)>=minimum)&&all(value(:)<=maximum) ...
    &&(isscalar(value)||numel(value)==nj),'JointOpeningProfile:JointParameter', ...
    'Invalid %s inventory or value.',name);
value=double(value(:)); if isscalar(value), value=repmat(value,nj,1); end
end
