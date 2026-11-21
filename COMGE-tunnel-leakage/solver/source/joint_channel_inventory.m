function [E,dEdT,info]=joint_channel_inventory(T,freeze,storage,c,L,rho,latent,cw,ci,rhoCcc,g)
%JOINT_CHANNEL_INVENTORY Equilibrium microchannel water/ice enthalpy.
% T is nj-by-nq. E is excess energy over displaced concrete per sample.
% It uses the r^2 VOLUME-weighted truncation, distinct from r^4 FLOW weight.
% Equal water/ice densities and fixed total saturated pore volume are used,
% consistently with the surrounding TH model. No permanent damage is added.
% A locally homogenized capillary network has, at reference temperature,
% Cref=rho*g*A*rQ^2*exp(-2*s^2)/(8*muRef*L*tortuosity).
% A is channel flow area per metre of tunnel axis [m2/m], not a porosity.
% Tortuosity is the geometric path-length ratio Lpath/L >= 1; consequently
% channel volume per metre of tunnel axis is A*Lpath [m3/m = m2].
% Thus specified conductance/radius/tortuosity define equivalent connected
% area and channel inventory. Tortuosity and the equivalent geometry remain
% uncalibrated scenario inputs, not inferred physical gasket measurements.
% Optional G defaults to 9.81 m/s2. L may be scalar or one value per joint.
% If ln(r)~N(m,s^2), rQ=exp(m+4s^2) while rV=exp(m+2s^2).
% Therefore rV=rQ*exp(-2s^2), and the connected-volume fraction remains
% liquid below Tbulk according to Phi((log(rcrit/rQ)+2s^2)/s).
assert(nargin==10||nargin==11,'joint_channel_inventory:InvalidInputCount', ...
    'Supply ten inputs and, optionally, gravity g.');
if nargin<11, g=9.81; end
assert(isnumeric(T)&&isreal(T)&&ismatrix(T)&&all(isfinite(T(:))), ...
    'joint_channel_inventory:InvalidTemperature', ...
    'T must be a finite real joint-by-sample matrix.');
T=double(T); [nj,nq]=size(T);
assert(nq>0||nj==0,'joint_channel_inventory:EmptyPath', ...
    'Every joint requires at least one temperature sample.');
freeze=joint_freeze_defaults(freeze);
assert(isstruct(storage)&&isscalar(storage)&& ...
    all(isfield(storage,{'enabled','tortuosity'})), ...
    'joint_channel_inventory:InvalidStorage', ...
    'Storage must supply enabled and geometric tortuosity.');
validateSwitch(storage.enabled,'storage.enabled');
validateSwitch(freeze.enabled,'freeze.enabled');
validateSwitch(freeze.iceEnabled,'freeze.iceEnabled');
c=expand(c,nj,false,'Conductance');
rq=expand(freeze.flowWeightedRadius_m,nj,true,'Radius');
[rqProfile,~]=joint_radius_profile(freeze,nj,nq);
s=expand(freeze.logRadiusStd,nj,true,'LogRadiusStd');
tau=expand(storage.tortuosity,nj,true,'Tortuosity');
L=expand(L,nj,true,'PathThickness');
assert(all(tau>=1),'joint_channel_inventory:InvalidTortuosity', ...
    'Geometric channel path-length ratio must be >=1.');
validatePositiveScalar(rho,'Density'); validatePositiveScalar(g,'Gravity');
validatePositiveScalar(latent,'LatentHeat'); validatePositiveScalar(cw,'WaterHeatCapacity');
validatePositiveScalar(ci,'IceHeatCapacity'); validatePositiveScalar(rhoCcc,'ConcreteHeatCapacity');
validatePositiveScalar(freeze.interfacialEnergy_J_m2,'InterfacialEnergy');
validatePositiveScalar(freeze.iceDensity_kg_m3,'IceDensity');
validatePositiveScalar(freeze.latentHeat_J_kg,'FreezingLatentHeat');
validatePositiveScalar(freeze.meltingTemperature_K,'MeltingTemperature');
assert(isnumeric(freeze.bulkFreezingTemperature_C)&& ...
    isreal(freeze.bulkFreezingTemperature_C)&& ...
    isscalar(freeze.bulkFreezingTemperature_C)&& ...
    isfinite(freeze.bulkFreezingTemperature_C), ...
    'joint_channel_inventory:InvalidBulkFreezingTemperature', ...
    'Bulk freezing temperature must be a finite real scalar.');
mu=water_viscosity_iapws(freeze.referenceTemperature_C);
area=zeros(nj,1); positive=c>0;
logArea=log(8)+log(mu)+log(L(positive))+log(c(positive)) ...
    +2*s(positive).^2+log(tau(positive))-log(rho)-log(g)-2*log(rq(positive));
area(positive)=exp(logArea);
pathLength=L.*tau; volume=area.*pathLength;
assert(all(isfinite(area))&&all(isfinite(volume)), ...
    'joint_channel_inventory:NonFiniteInventory', ...
    'Equivalent channel area or volume overflows; check physical inputs.');
fv=ones(nj,nq); df=zeros(nj,nq);
if freeze.enabled&&freeze.iceEnabled
    depression=freeze.bulkFreezingTemperature_C-T;
    cold=depression>0;
    logBeta=log(2)+log(freeze.interfacialEnergy_J_m2)+ ...
        log(freeze.meltingTemperature_K)-log(freeze.iceDensity_kg_m3)- ...
        log(freeze.latentHeat_J_kg);
    rqGrid=rqProfile; sGrid=repmat(s,1,nq);
    z=(logBeta-log(depression(cold))-log(rqGrid(cold))+ ...
        2*sGrid(cold).^2)./sGrid(cold);
    fv(cold)=.5*erfc(-z/sqrt(2));
    df(cold)=exp(-.5*z.^2-.5*log(2*pi)-log(sGrid(cold))-log(depression(cold)));
end
sampleVolume=repmat(volume/nq,1,nq);
E=zeros(nj,nq); dEdT=zeros(nj,nq); active=positive&logical(storage.enabled);
E(active,:)=sampleVolume(active,:).*(rho*(latent*fv(active,:)+ ...
    (ci+(cw-ci)*fv(active,:)).*T(active,:))-rhoCcc*T(active,:));
dEdT(active,:)=sampleVolume(active,:).*(rho*(ci+(cw-ci)*fv(active,:)+ ...
    (latent+(cw-ci)*T(active,:)).*df(active,:))-rhoCcc);
assert(all(isfinite(E(:)))&&all(isfinite(dEdT(:))), ...
    'joint_channel_inventory:NonFiniteEnergy', ...
    'Channel energy or its temperature derivative is not finite.');
if ~storage.enabled
    E(:)=0; dEdT(:)=0; volume(:)=0; area(:)=0;
end
info=struct('connectedArea_m2_per_m',area,'waterVolume_m2',volume, ...
    'liquidVolumeFraction',fv,'temperatureDerivative_per_C',df, ...
    'pathLength_m',pathLength,'volumeWeightedRadius_m',rq.*exp(-2*s.^2));
info.localVolumeWeightedRadius_m=rqProfile.*exp(-2*s.^2);
end

function x=expand(x,n,positive,name)
assert(isnumeric(x)&&isreal(x)&&(isvector(x)||isempty(x)) ...
    &&all(isfinite(x(:)))&&(isscalar(x)||numel(x)==n), ...
    ['joint_channel_inventory:Invalid' name], ...
    '%s must be finite, real, and scalar or one value per joint.',name);
x=double(x(:));
if positive
    assert(all(x>0),['joint_channel_inventory:Invalid' name], ...
        '%s must be strictly positive, independently of conductance.',name);
else
    assert(all(x>=0),['joint_channel_inventory:Invalid' name], ...
        '%s must be nonnegative.',name);
end
if isscalar(x),x=repmat(x,n,1);end
end

function validatePositiveScalar(x,name)
assert(isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0, ...
    ['joint_channel_inventory:Invalid' name], ...
    '%s must be a positive finite real scalar.',name);
end

function validateSwitch(x,name)
assert((isnumeric(x)||islogical(x))&&isreal(x)&&isscalar(x)&& ...
    isfinite(x)&&any(x==[0 1]),'joint_channel_inventory:InvalidSwitch', ...
    '%s must be scalar logical or numeric 0/1.',name);
end
