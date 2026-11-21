function [mu,dmu]=water_viscosity_iapws(T)
%WATER_VISCOSITY_IAPWS Liquid-water viscosity at about 0.1 MPa, Pa s.
% IAPWS SR6-08(2011), Eq.7/Table5; -20..110 C including supercooled liquid.
% https://iapws.org/technical-guidance/release/LiquidWater
assert(isnumeric(T)&&isreal(T)&&all(isfinite(T(:)))&& ...
    all(T(:)>=-20&T(:)<=110),'water_viscosity_iapws:range', ...
    'Liquid viscosity correlation requires -20 <= T <= 110 C.');
tk=T+273.15; mu=zeros(size(T));
wantDerivative=nargout>1;
if wantDerivative, dmu=mu; end
a=[280.68,511.45,61.131,.45903]; b=[-1.9,-7.7,-19.6,-40];
% Share log(TK/300) across the four powers of the same positive base.
% This is the same correlation, with differences at roundoff level.
logRatio=log(tk/300);
for k=1:4
    term=1e-6*a(k)*exp(b(k)*logRatio);
    mu=mu+term;
    if wantDerivative, dmu=dmu+term*b(k)./tk; end
end
end
