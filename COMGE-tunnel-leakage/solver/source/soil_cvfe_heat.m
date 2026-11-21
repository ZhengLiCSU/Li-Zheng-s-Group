function [rHeat,DT,DF,diagnostic] = soil_cvfe_heat(T,F,G,rho,L,cw)
%SOIL_CVFE_HEAT Conservative upstream enthalpy transport on CVFE subfaces.
% [rHeat,DT,DF,diagnostic]=soil_cvfe_heat(Tglobal,F,G,rho,L,cw) returns the
% GLOBAL temperature residual, partial dResidual/dT holding F fixed, and
% partial dResidual/dF. Interior subface contributions sum exactly to zero.
% Liquid specific enthalpy uses the core reference rho*(L+cw*Tup). Heat is
% selected by the actual Darcy normal flow, never by a cotangent edge sign.
%
% This is first-order nodal upstream transport on the barycentric dual;
% it is not a higher-order reconstruction or a mixed RT0 solution. Boundary
% heat fluxes and separate joint-channel inventories are caller terms.
% For nonlinear conductivity the full tangent is DT+DF*dF/dT, DF*dF/dh.
% At F=0 the from-node donor gives a selected one-sided derivative; the
% residual itself is continuous. A finite-difference test across that switch
% cannot require a unique classical derivative for unequal donor T.

validateattributes(T,{'numeric'},{'vector','real','finite','numel',G.N}, ...
    mfilename,'T',1);
validateattributes(F,{'numeric'},{'vector','real','finite','numel',G.nf}, ...
    mfilename,'F',2);
validateattributes(rho,{'numeric'},{'scalar','real','finite','positive'}, ...
    mfilename,'rho',4);
validateattributes(L,{'numeric'},{'scalar','real','finite'},mfilename,'L',5);
validateattributes(cw,{'numeric'},{'scalar','real','finite','nonnegative'}, ...
    mfilename,'cw',6);
T=T(:); F=F(:); donor=G.toGlobal;
forward=F>=0; donor(forward)=G.fromGlobal(forward);
enthalpy=rho*(L+cw*T(donor)); heatFlow=enthalpy.*F;
globalInc=G.soilSelect'*G.inc;
rHeat=globalInc*heatFlow;
if nargout>=2
    donorSlope=sparse((1:G.nf)',donor,rho*cw*F,G.nf,G.N);
    DT=globalInc*donorSlope;
end
if nargout>=3
    DF=globalInc*spdiags(enthalpy,0,G.nf,G.nf);
end
if nargout>=4
    diagnostic=struct('upstreamGlobalNode',donor, ...
        'upstreamTemperature',T(donor),'liquidEnthalpy',enthalpy, ...
        'advectiveEnthalpyFlux',heatFlow);
end
end
