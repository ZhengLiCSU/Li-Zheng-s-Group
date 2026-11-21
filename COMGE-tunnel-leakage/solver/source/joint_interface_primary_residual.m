function [residual,tangent,flow,state,x] = joint_interface_primary_residual(z,old,temperature,dt,settings,p)
%JOINT_INTERFACE_PRIMARY_RESIDUAL Pressure coordinates for the same equations.
% z = [gap increment m, pressure Pa, ice MASS fraction]. Mass is reconstructed
% from the exact linear-compressibility volume relation; the mass PDE remains
% in the residual. This avoids pressure cancellation for long time steps.
arguments
    z (:,3) double {mustBeReal,mustBeFinite}
    old (:,4) double {mustBeReal,mustBeFinite}
    temperature (1,1) double {mustBePositive,mustBeFinite}
    dt (1,1) double {mustBePositive,mustBeFinite}
    settings (1,1) struct
    p (1,1) struct
end
n = size(z,1); fraction = z(:,3); pressure = z(:,2);
vw = 1/p.rhoWater; vi = 1/p.rhoIce;
cw = vw/p.bulkWater_Pa; ci = vi/p.bulkIce_Pa;
v = vw+fraction*(vi-vw); c = cw+fraction*(ci-cw);
denominator = v-c.*pressure;
mass = (p.aperture0_m+p.wetAreaFraction*z(:,1))./denominator;
x = [z(:,1),mass,fraction.*mass];
[residual,original,flow,state] = joint_interface_residual(x,old,temperature,dt,settings,p,pressure);
massGap = p.wetAreaFraction./denominator;
massPressure = mass.*c./denominator;
massFraction = -mass.*((vi-vw)-(ci-cw)*pressure)./denominator;
transform = sparse(3*n,3*n);
for k = 1:n
    ids = 3*k-2:3*k;
    transform(ids,ids) = [1,0,0;massGap(k),massPressure(k),massFraction(k); ...
        fraction(k)*massGap(k),fraction(k)*massPressure(k),mass(k)+fraction(k)*massFraction(k)];
end
tangent = original*transform;
end
