function v = premelting_film_potential(h,f)
%PREMELTING_FILM_POTENTIAL One wall-film-ice potential and exact derivatives.
% W=C exp(-k h)-A/(12 pi h^2); Pi=-W_h. Positive A is ATTRACTIVE.
% Stability requires W_hh>0. A positive Hamaker constant must never be
% silently treated as a repulsive A/(6*pi*h^3) pressure in this convention.
arguments
    h (:,1) double {mustBeFinite,mustBePositive}
    f (1,1) struct = premelting_film_parameters()
end
C=f.shortRangeEnergy_Jm2; k=f.inverseLength_minv; A=f.hamaker_J;
assert(isfinite(C)&&C>=0&&isfinite(k)&&k>0&&isfinite(A), ...
    'Premelting:Material','Invalid interface potential.');
e=C*exp(-k*h);
W=e-A./(12*pi*h.^2);
Pi=k*e-A./(6*pi*h.^3);
PiH=-k^2*e+A./(2*pi*h.^4);
v=struct('energy_Jm2',W,'pressure_Pa',Pi,'pressureDerivative_Pam',PiH, ...
    'energyDerivative_Jm3',-Pi,'curvature_Jm4',-PiH);
end
