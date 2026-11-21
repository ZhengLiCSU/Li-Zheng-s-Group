function [mu,dT,dRho] = water_viscosity_iapws2008(T,rho)
%WATER_VISCOSITY_IAPWS2008 Bulk-water viscosity and exact analytic derivatives.
% IAPWS R12-08, Eq.11-12/Table1-2, critical enhancement omitted (Sec.2.8).
% https://iapws.org/documents/release/viscosity.download
% The bulk correlation is not a measurement of nanofilm viscosity. Negative
% liquid pressure and molecular confinement are model extrapolations.
arguments
    T (:,1) double {mustBeFinite,mustBePositive}
    rho (:,1) double {mustBeFinite,mustBeNonnegative}
end
assert(numel(T)==numel(rho) && all(T>=250 & T<=1173.15), ...
    'WaterViscosity:Domain','Supported temperature range is 250..1173.15 K with matched input sizes.');
assert(~any(T>645.91 & T<650.77 & rho>245.8 & rho<405.3), ...
    'WaterViscosity:Domain','Near-critical enhancement is not implemented.');
h=[1.67752,2.20462,0.6366564,-0.241605];
H=[.520094,.222531,-.281378,.161913,-.0325372,0,0; ...
   .0850895,.999115,-.906851,.257399,0,0,0; ...
   -1.08374,1.88797,-.772479,0,0,0,0; ...
   -.289555,1.26613,-.489837,0,.0698452,0,-.00435673; ...
   0,0,-.257040,0,0,.00872102,0; ...
   0,.120573,0,0,0,0,-.000593264];
tr=T/647.096; rr=rho/322; x=1./tr-1; y=rr-1;
den=zeros(size(T)); dden=den;
for i=0:3
    den=den+h(i+1)*tr.^(-i);
    if i>0, dden=dden-i*h(i+1)*tr.^(-i-1)/647.096; end
end
S=zeros(size(T)); Sx=S; Sy=S;
for i=0:5
    for j=0:6
        a=H(i+1,j+1); if a==0, continue; end
        S=S+a*x.^i.*y.^j;
        if i>0, Sx=Sx+i*a*x.^(i-1).*y.^j; end
        if j>0, Sy=Sy+j*a*x.^i.*y.^(j-1); end
    end
end
mu=1e-6*100*sqrt(tr)./den.*exp(rr.*S);
dT=mu.*(0.5./T-dden./den-rr.*Sx.*647.096./T.^2);
dRho=mu.*(S+rr.*Sy)/322;
assert(all(isfinite(mu)&mu>0),'WaterViscosity:Domain','Correlation produced invalid viscosity.');
end
