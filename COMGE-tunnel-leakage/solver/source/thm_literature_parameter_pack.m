function [model,evidence]=thm_literature_parameter_pack(name)
%THM_LITERATURE_PARAMETER_PACK Traceable SOIL overrides, not site calibration.
% name='kanagawa_th' (default) or 'silt_retention'. Values are read from
% Dall'Amico et al., TC5(2011)469, Section7.2 and Table1, respectively.
% Returns minimal production overrides and a separate provenance ledger.
% No rock fatigue or concrete-crack parameters are assigned to EPDM joints.
% The tunnel boundary/mechanics/channel case remains explicitly assumed.
if nargin<1,name='kanagawa_th';end
assert((ischar(name)&&isrow(name))||(isstring(name)&&isscalar(name)), ...
    'THMLiterature:Name','Supply one named literature pack.');
source="https://doi.org/10.5194/tc-5-469-2011";
switch char(name)
    case 'kanagawa_th'
        % mm^-1 -> m^-1; mm/s -> m/s. Solid volumetric capacity is rhoS*cs.
        model=struct('porosity',.535,'thetaR',.05,'vgAlpha',1.11, ...
            'vgN',1.48,'Ksat',3.2e-6,'ks',2.5,'rhoS',2650,'cs',2.3e6/2650);
        keys={'porosity';'thetaR';'vgAlpha';'vgN';'Ksat';'ks';'rhoS*cs'};
        values=[.535;.05;1.11;1.48;3.2e-6;2.5;2.3e6];
        units={'1';'1';'m^-1';'1';'m/s';'W/(m K)';'J/(m^3 K)'};
        original=[.535;.05;1.11e-3;1.48;.0032;2.5;2.3e6];
        originalUnits={'1';'1';'mm^-1';'1';'mm/s';'W/(m K)';'J/(m^3 K)'};
        location='Section7.2, Kanagawa sandy loam benchmark';
        conversion={'identity';'identity';'multiply1000';'identity'; ...
            'multiply0.001';'identity';'rhoS*cs=published capacity; rhoS split assumed'};
        missing={'Soil mechanical parameters; saturated storage; reference temperature of Ksat', ...
            'Joint conductance/radius distribution; gasket failure; cyclic hydraulic deterioration', ...
            'Tunnel air/groundwater histories; tunnel and soil geometry'};
    case 'silt_retention'
        % This four-parameter set is NOT a complete experimental soil.
        model=struct('porosity',.49,'thetaR',.05,'vgAlpha',.65,'vgN',1.67);
        keys={'porosity';'thetaR';'vgAlpha';'vgN'};values=[.49;.05;.65;1.67];
        units={'1';'1';'m^-1';'1'};original=[.49;.05;6.5e-4;1.67];
        originalUnits={'1';'1';'mm^-1';'1'};location='Table1, silt retention parameters';
        conversion={'identity';'identity';'multiply1000';'identity'};
        missing={'Ksat and thermal properties are NOT specified by this retention table', ...
            'All tunnel, joint, cycle-damage and mechanical inputs remain scenario choices'};
    otherwise
        error('THMLiterature:Name','Unknown pack: %s.',char(name));
end
ledger=table(string(keys),values,string(units),original,string(originalUnits), ...
    string(conversion),repmat(source,numel(keys),1), ...
    repmat(string(location),numel(keys),1), ...
    'VariableNames',{'parameter','modelValue','unit','publishedValue', ...
    'publishedUnit','conversion','source','location'});
evidence=struct('name',char(name),'schemaVersion',1,'ledger',ledger, ...
    'overrides',model,'missingCalibration',{missing}, ...
    'scope','Literature soil inputs mapped to present constitutive model; no reproduction of full column experiment or tunnel calibration.', ...
    'isSiteCalibration',false,'isJointDamageCalibration',false);
end
