function f = premelting_film_parameters(substrate)
%PREMELTING_FILM_PARAMETERS Published model-wall coefficients, NOT EPDM.
% Baran, Llombart & MacDowell, J Phys Chem C (2025),
% doi:10.1021/acs.jpcc.4c07328; public preprint arXiv:2410.21893v1,
% Table 2, wall rows f=1..4. Coefficients are frozen at their tabulated state.
% The paper's TIP4P/Ice melting point is 270 K, not real water's triple point.
% Use to reproduce the INTERFACE POTENTIAL, not its entire MD experiment.
arguments
    substrate (1,1) double {mustBeInteger,mustBePositive} = 4
end
assert(substrate<=4,'Premelting:Material','Only the four published model walls are provided.');
C=[31.330,30.322,66.079,40.969]*1e-3;
kappa=[0.323,0.258,0.289,0.166]*1e10;
A=[2.27,1.36,0.66,0.071]*1e-21;
he=[31,44,45,100]*1e-10;
f=struct('shortRangeEnergy_Jm2',C(substrate),'inverseLength_minv',kappa(substrate), ...
    'hamaker_J',A(substrate),'minimumFilm_m',1e-9,'maximumFilm_m',1e-6, ...
    'maximumFilmFraction',0.1,'tabulatedCoexistenceFilm_m',he(substrate), ...
    'modelWallIndex',substrate,'source','10.1021/acs.jpcc.4c07328 / arXiv:2410.21893v1 Table 2', ...
    'scope','Illustrative frozen-coefficient model-wall potential; not calibrated concrete or gasket.');
end
