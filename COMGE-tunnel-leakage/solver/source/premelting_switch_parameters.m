function f = premelting_switch_parameters()
%PREMELTING_SWITCH_PARAMETERS Illustrative three-domain slit closure.
% Published W(h) coefficients are inherited from model wall f=4. The added
% constant surface energies are illustrative, NOT measured EPDM parameters.
% Parallel liquid / wet-ice / dry-vapour columns neglect lateral interfaces.
f=premelting_film_parameters(4);
f.iceWaterSurfaceEnergy_Jm2=0.030;
f.wallVaporExcess_Jm2=0.0756*cosd(50);
f.scope='Equilibrium parallel-domain convexification; no nucleation kinetics, air or lateral interfaces.';
f.spinodalFilm_m=fzero(@(h)curvature(h,f),[f.minimumFilm_m,f.maximumFilm_m],optimset('TolX',1e-22));
end

function c=curvature(h,f)
v=premelting_film_potential(h,f); c=v.curvature_Jm4;
end
