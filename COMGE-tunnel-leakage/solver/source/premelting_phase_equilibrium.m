function e = premelting_phase_equilibrium(pl,T,p,f)
%PREMELTING_PHASE_EQUILIBRIUM Distinct bulk chemical potentials + stable film.
% Existing ice, pure water, flat interfaces, no solute, no nucleation or
% vapour. Reject unsupported states; never clamp pressure or film thickness.
% Inputs pl gauge Pa, T kelvin. Columns of derivatives: [pl,T].
arguments
    pl (:,1) double {mustBeFinite}
    T (:,1) double {mustBeFinite}
    p (1,1) struct = joint_threephase_parameters()
    f (1,1) struct = premelting_film_parameters()
end
assert(numel(pl)==numel(T),'Premelting:Shape','Pressure and temperature sizes differ.');
ql=joint_threephase_gibbs(pl,T,p);
assert(all(ql.g(:,3)>=ql.g(:,1)),'Premelting:Cavitation', ...
    'The wet-ice branch cannot represent stable vapour; use a coupled three-phase formulation.');
Tm=p.meltingTemperature_K;
base=p.heatIce_JkgK*(T-Tm-T.*log(T/Tm))+ ...
    p.latentHeat_Jkg*(T/Tm-1)+p.iceReferenceShift_Jkg;
vi0=1/p.rhoIce; ci=vi0/p.bulkIce_Pa;
D=ql.g(:,1)-base; disc=vi0^2-2*ci*D;
assert(all(disc>0),'Premelting:EOS','Ice EOS has left its small-compression branch.');
% Cancellation-free root of vi0*pi - ci*pi^2/2 = D.
piIce=2*D./(vi0+sqrt(disc));
split=piIce-pl;
assert(all(split>0),'Premelting:IceAbsent', ...
    'A positive ice-liquid pressure split is required for this pre-existing wet-ice branch.');
qi=joint_threephase_gibbs(piIce,T,p);
lo=f.minimumFilm_m; hi=f.maximumFilm_m;
vlo=premelting_film_potential(lo,f); vhi=premelting_film_potential(hi,f);
assert(vlo.pressureDerivative_Pam<0 && vhi.pressure_Pa<min(split), ...
    'Premelting:Bracket','The requested stable film is not bracketed by the documented domain.');
assert(all(split<vlo.pressure_Pa),'Premelting:MolecularLimit', ...
    'Required film is thinner than the explicit continuum limit; no film floor is applied.');
h=zeros(size(pl));
for k=1:numel(h)
    h(k)=fzero(@(x)filmEquation(x,split(k),f),[lo,hi],optimset('TolX',1e-22));
end
v=premelting_film_potential(h,f);
assert(all(v.curvature_Jm4>0),'Premelting:UnstableFilm','Stationary film is unstable.');
dpi=[ql.volume(:,1)./qi.volume(:,2), ...
    (qi.entropy(:,2)-ql.entropy(:,1))./qi.volume(:,2)];
dh=(dpi-[ones(size(pl)),zeros(size(pl))])./v.pressureDerivative_Pam;
e=struct('liquidPressure_Pa',pl,'icePressure_Pa',piIce,'wallPressure_Pa',piIce, ...
    'temperature_K',T,'film_m',h,'pressureSplit_Pa',split,'potential',v, ...
    'liquid',ql,'ice',qi,'icePressureDerivative',dpi,'filmDerivative',dh, ...
    'chemicalPotential_Jkg',ql.g(:,1),'affinity_Jkg',qi.g(:,2)-ql.g(:,1));
end

function r=filmEquation(h,split,f)
v=premelting_film_potential(h,f); r=v.pressure_Pa-split;
end
