function [s,d] = thermo_interface_state(z,model,p)
%THERMO_INTERFACE_STATE Exact static condensation of LINEAR host elasticity.
% z=[pressure;ice mass fraction;joint temperature;solid temperature]. The
% host includes supported Q4 elasticity and reversible normal restraint.
% Thermal expansion and pressure are coupled to the SAME gap operator.
n = numel(model.weight); ns = numel(model.capacity); nz = 3*n+ns;
assert(numel(z)==nz,'ThermoInterface:Shape','State size must match host and interface.');
ip=1:n; iff=n+(1:n); it=2*n+(1:n); its=3*n+(1:ns);
pressure=z(ip); fraction=z(iff); temperature=z(it); solidT=z(its);
assert(all(solidT>0),'JointNonisothermal:Temperature','Solid absolute temperatures must be positive.');
u = model.displacementPressure*pressure+model.displacementTemperature*(solidT-model.referenceTemperature_K);
gap = model.gapOperator*u;
vw=1/p.rhoWater; vi=1/p.rhoIce; cw=vw/p.bulkWater_Pa; ci=vi/p.bulkIce_Pa;
c=cw+fraction*(ci-cw); v=vw+fraction*(vi-vw); denominator=v-c.*pressure;
mass=(p.aperture0_m+p.wetAreaFraction*gap)./denominator; ice=fraction.*mass;
pf=p; pf.normalStiffness_Pam=0; pf.contactStiffness_Pam=0; pf.externalTraction_Pa=0;
s=joint_helmholtz_state([gap,mass,ice,temperature],pf,true,pressure);
s.x=[gap,mass,ice,temperature]; s.u=u; s.solidTemperature_K=solidT;
s.pressure=pressure; s.fraction=fraction;
if nargout<2, return; end
du=zeros(size(model.stiffness,1),nz);
du(:,ip)=model.displacementPressure; du(:,its)=model.displacementTemperature;
dg=model.gapOperator*du;
dp=zeros(n,nz); dp(:,ip)=eye(n); df=zeros(n,nz); df(:,iff)=eye(n);
dT=zeros(n,nz); dT(:,it)=eye(n);
dm=(p.wetAreaFraction./denominator).*dg+(mass.*c./denominator).*dp- ...
    (mass.*((vi-vw)-(ci-cw)*pressure)./denominator).*df;
di=fraction.*dm+mass.*df;
dl=(vw-cw*pressure).*(dm-di)-(mass-ice)*cw.*dp;
dv=p.wetAreaFraction*dg;
de=s.energyGradient(:,1).*dg+s.energyGradient(:,2).*dm+ ...
    s.energyGradient(:,3).*di+s.energyGradient(:,4).*dT;
dA=zeros(n,nz);
for k=1:n
    dA(k,:)=s.hessian(3,:,k)*[dg(k,:);dm(k,:);di(k,:);dT(k,:)];
end
d=struct('u',du,'gap',dg,'pressure',dp,'fraction',df,'temperature',dT, ...
    'mass',dm,'ice',di,'liquidVolume',dl,'iceVolume',dv-dl,'energy',de,'affinity',dA);
end
