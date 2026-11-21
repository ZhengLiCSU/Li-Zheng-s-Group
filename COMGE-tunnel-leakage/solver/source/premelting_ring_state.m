function [s,d] = premelting_ring_state(pl,T,Ts,model,p,f,history)
%PREMELTING_RING_STATE Two-pressure film thermodynamics coupled to Q4 ring.
% An equilibrium constitutive bridge, NOT a transient flow/heat solution.
% Ice pressure acts on the wall; liquid pressure drives water transport.
e=premelting_phase_equilibrium(pl,T,p,f);
m=solve_cohesive_ring_mechanics(e.wallPressure_Pa,Ts,model,p,history);
gap=model.gapOperator*m.u; b=p.aperture0_m+p.wetAreaFraction*gap;
[q,v]=premelting_slit_state(pl,T,b,p,f);
res=model.stiffness*m.u+model.gapOperator'*(model.weight.*m.local.traction_Pa)- ...
    model.thermalLoad*(Ts-model.referenceTemperature_K)- ...
    model.gapOperator'*(model.weight*p.wetAreaFraction.*e.wallPressure_Pa);
assert(norm(res,inf)<1e-4,'Premelting:Mechanics','Ice wall-load balance failed.');
strain=reshape(model.strainOperator*m.u,3,[])';
meanE=(strain(:,1)+strain(:,2))/2; radius=hypot((strain(:,1)-strain(:,2))/2,strain(:,3)/2);
maximumStrain=max(abs([meanE+radius;meanE-radius]));
assert(maximumStrain<0.05,'Premelting:Strain','Small strain domain exceeded.');
s=struct('slit',q,'u',m.u,'gap_m',gap,'cohesive',m.local, ...
    'solidTemperature_K',Ts,'mechanicalResidual_Nm',res,'maximumPrincipalStrain',maximumStrain);
if nargout<2, return; end
up=m.pressureDerivative.*e.icePressureDerivative(:,1)';
ut=m.pressureDerivative.*e.icePressureDerivative(:,2)';
uTs=m.temperatureDerivative;
db=p.wetAreaFraction*model.gapOperator*[up,ut,uTs]; n=numel(pl); ns=numel(Ts);
dM=[diag(v.mass(:,1)),diag(v.mass(:,2)),zeros(n,ns)]+v.mass(:,3).*db;
d=struct('u',[up,ut,uTs],'aperture',db,'mass',dM, ...
    'energy',[diag(v.energy(:,1)),diag(v.energy(:,2)),zeros(n,ns)]+v.energy(:,3).*db, ...
    'columns','[liquid pressure; joint temperature; solid temperature]');
end
