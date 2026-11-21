function bc = thermo_interface_boundary(model,outerTemperature,innerTemperature,pressureDrop,sealed)
%THERMO_INTERFACE_BOUNDARY Explicit thermal baths and liquid reservoirs.
arguments
    model (1,1) struct
    outerTemperature (1,1) double {mustBePositive,mustBeFinite}
    innerTemperature (1,1) double {mustBePositive,mustBeFinite}
    pressureDrop (1,1) double {mustBeNonnegative,mustBeFinite} = 2e5
    sealed (1,1) logical = false
end
ns=numel(model.capacity); denominator=model.bathConductance;
Tb=repmat(model.referenceTemperature_K,ns,1); active=denominator>0;
Tb(active)=(model.innerBathConductance(active)*innerTemperature+ ...
    model.outerBathConductance(active)*outerTemperature)./denominator(active);
bc=struct('reservoirPressure_Pa',repmat([pressureDrop;0],8,1), ...
    'reservoirTemperature_K',repmat([outerTemperature;innerTemperature],8,1), ...
    'sealedMass',repmat(sealed,16,1),'insulatedHeat',false(16,1), ...
    'outletSeepage',repmat([false;true],8,1),'solidBathTemperature_K',Tb);
end
