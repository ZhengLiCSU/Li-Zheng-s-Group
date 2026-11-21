function [s,d] = premelting_switch_ring_state(z,model,p,f,history,loadingAtHistory)
%PREMELTING_SWITCH_RING_STATE Independent liquid chemical pressure and wall load.
% z=[pl/1e6; theta_l; theta_i; theta_v; Tjoint; Pwall/1e6; Tsolid].
% Mechanics sees the conjugate WALL pressure; liquid transport sees pl.
if nargin<6, loadingAtHistory=[]; end
model=prepare_premelting_solver(model,p);
n=numel(model.weight); ns=numel(model.capacity); nz=6*n+ns;
assert(numel(z)==nz,'PhaseSwitch:Shape','Incorrect state size.');
pl=z(1:n)*1e6; theta=reshape(z(n+(1:3*n)),n,3);
T=z(4*n+(1:n)); wall=z(5*n+(1:n))*1e6; Ts=z(6*n+(1:ns));
assert(all(Ts>=240 & Ts<=310),'PhaseSwitch:Domain','Solid temperature domain exceeded.');
mech=solve_cohesive_ring_mechanics(wall,Ts,model,p,history,loadingAtHistory,nargout>1);
gap=model.gapOperator*mech.u; b=p.aperture0_m+gap;
% Invalid Newton trials must return the domain identifier recognized by
% the line search, before MATLAB argument validators abort the solver.
assert(all(isfinite(b)&b>2*f.minimumFilm_m),'PhaseSwitch:Domain','Trial aperture is outside the molecular-film domain.');
if nargout>1
    [local,ld]=premelting_switch_state(pl,T,b,theta,p,f);
else
    local=premelting_switch_state(pl,T,b,theta,p,f);
end
s=local; s.local=local; s.u=mech.u; s.solidTemperature_K=Ts;
s.x=[gap,local.mass_kgm2,local.phaseMass_kgm2(:,2),T];
s.pressure=wall; s.liquidPressure=pl; s.cohesive=mech.local;
s.liquidVolume_m=local.phaseVolume_m(:,1); s.iceVolume_m=local.phaseVolume_m(:,2); s.vaporVolume_m=local.phaseVolume_m(:,3);
if nargout<2, return; end
dp=zeros(n,nz); dp(:,1:n)=1e6*eye(n);
dT=zeros(n,nz); dT(:,4*n+(1:n))=eye(n);
dw=zeros(n,nz); dw(:,5*n+(1:n))=1e6*eye(n);
% The pressure map is supported on wall-pressure columns only. Avoid a
% dense multiplication by an otherwise all-zero derivative matrix.
du=zeros(numel(mech.u),nz);
du(:,5*n+(1:n))=1e6*mech.pressureDerivative;
du(:,6*n+(1:ns))=mech.temperatureDerivative;
dg=model.gapOperator*du;
areas=zeros(n,3,nz);
for k=1:3
    a=zeros(n,nz); a(:,k*n+(1:n))=eye(n);
    areas(:,k,:)=a;
end
d=struct('u',du,'gap',dg,'pressure',dw,'liquidPressure',dp,'temperature',dT,'areas',areas);
% The six local arguments have disjoint global supports: local liquid
% pressure, local joint T, mechanical gap, and three local area fractions.
% Bind only their nonzero blocks instead of multiplying six zero-filled
% arrays for every property. All entries are the same analytic chain rule.
mapping=struct('n',n,'nz',nz,'gap',dg(:,5*n+1:end));
names={'mass','energy','entropy','freeEnergy','wallPressure','chemicalMobility','volumeConductance'};
for k=1:numel(names)
    d.(names{k})=bindLocalDerivative(ld.(names{k}),mapping);
end
d.grandPotential=zeros(n,3,nz); volumes=zeros(n,3,nz);
for k=1:3
    v=bindLocalDerivative(reshape(ld.phaseVolume(:,k,:),n,6),mapping);
    a=bindLocalDerivative(reshape(ld.grandPotential(:,k,:),n,6),mapping);
    volumes(:,k,:)=v; d.grandPotential(:,k,:)=a;
end
d.liquidVolume=reshape(volumes(:,1,:),n,nz); d.iceVolume=reshape(volumes(:,2,:),n,nz); d.vaporVolume=reshape(volumes(:,3,:),n,nz);
end

function value=bindLocalDerivative(local,mapping)
n=mapping.n; nz=mapping.nz; row=(1:n)';
value=zeros(n,nz);
value(sub2ind([n,nz],row,row))=1e6*local(:,1);
value(sub2ind([n,nz],row,4*n+row))=local(:,2);
for phase=1:3
    value(sub2ind([n,nz],row,phase*n+row))=local(:,3+phase);
end
value(:,5*n+1:end)=local(:,3).*mapping.gap;
end
