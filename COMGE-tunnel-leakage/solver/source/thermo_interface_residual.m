function [r,J,a] = thermo_interface_residual(z,old,dt,model,bc,p)
%THERMO_INTERFACE_RESIDUAL Nonisothermal interface + thermoelastic host.
% Independent joint T, pressure and ice fraction. Host temperatures are
% unknowns. Midpoint pressure work and reciprocal thermoelastic heating
% make TOTAL internal energy conservative; no reconstructed thermostat heat.
% Linear mechanical condensation is exact, not a staggered approximation.
[s,d]=thermo_interface_state(z,model,p); s0=thermo_interface_state(old,model,p);
n=numel(model.weight); ns=numel(model.capacity); nz=numel(z); w=model.weight;
ip=1:n; iff=n+(1:n); it=2*n+(1:n); its=3*n+(1:ns);
vw=1/p.rhoWater; cw=vw/p.bulkWater_Pa;
vl=vw-cw*s.pressure;
mob=s.liquidVolume_m.^3./(12*p.viscosity_Pas*vl.^2*p.wetAreaFraction^2);
dmob=s.liquidVolume_m.^2./(4*p.viscosity_Pas*vl.^2*p.wetAreaFraction^2).*d.liquidVolume+ ...
    s.liquidVolume_m.^3*cw./(6*p.viscosity_Pas*vl.^3*p.wetAreaFraction^2).*d.pressure;
kappa=model.conductivityWater*s.liquidVolume_m+model.conductivityIce*s.iceVolume_m;
dkappa=model.conductivityWater*d.liquidVolume+model.conductivityIce*d.iceVolume;
nf=size(model.faces,1); flux=zeros(nf,2); dmass=zeros(nf,nz); denergy=dmass;
entropyFace=zeros(nf,1); boundaryEntropy=0;
for f=1:nf
    l=model.faces(f,1); q=model.faces(f,2); distance=model.faceDistance(f);
    [left,right,dl,dr,m,dm,k,dk]=faceData(l,q,f,s,d,mob,dmob,kappa,dkappa,model,bc,distance,nz);
    [value,derivative,e]=joint_nonisothermal_face(left,right,m,k,p);
    flux(f,:)=value';
    chain=[dl;dr;dm;dk]; gradient=derivative*chain;
    dmass(f,:)=gradient(1,:); denergy(f,:)=gradient(2,:);
    entropyFace(f)=e.entropyProduction_WmK;
    if l==0 || q==0
        b=model.boundaryID(f); Tb=bc.reservoirTemperature_K(b); pb=bc.reservoirPressure_Pa(b);
        mu=p.heatWater_JkgK*(Tb-p.meltingTemperature_K-Tb*log(Tb/p.meltingTemperature_K))+vw*pb-0.5*cw*pb^2;
        boundaryEntropy=boundaryEntropy+(1-2*(q==0))*(value(2)-mu*value(1))/Tb;
    end
end
% Trace heat exchange is a positive node-to-cell port graph (lumped trace
% quadrature), preserving both equal/opposite energy and port entropy.
H=model.heatExchange; solidT=s.solidTemperature_K; T=s.x(:,4);
referenceT=model.referenceTemperature_K;
jointHeat=model.exchangeRowSum.*(T-referenceT)-H*(solidT-referenceT);
solidHeat=model.exchangeColumnSum.*(solidT-referenceT)-H'*(T-referenceT);
djointHeat=model.exchangeRowSum.*d.temperature; djointHeat(:,its)=djointHeat(:,its)-H;
dsolidHeat=-H'*d.temperature; dsolidHeat(:,its)=dsolidHeat(:,its)+diag(model.exchangeColumnSum);
deltaGap=s.x(:,1)-s0.x(:,1); pAverage=0.5*(s.pressure+s0.pressure);
pressureWork=p.wetAreaFraction*pAverage.*deltaGap;
dwork=p.wetAreaFraction*(0.5*deltaGap.*d.pressure+pAverage.*d.gap);
massResidual=w.*(s.x(:,2)-s0.x(:,2))+dt*model.divergence*flux(:,1);
energyResidual=w.*(s.internalEnergy_Jm2-s0.internalEnergy_Jm2+pressureWork)+ ...
    dt*(model.divergence*flux(:,2)+jointHeat);
massJac=w.*d.mass+dt*model.divergence*dmass;
energyJac=w.*(d.energy+dwork)+dt*(model.divergence*denergy+djointHeat);
% Reciprocal term T d(thermalLoadMap' U) is REQUIRED by the solid free energy.
deltaU=s.u-s0.u; strainEntropy=model.thermalLoad'*deltaU;
averageSolidT=0.5*(solidT+s0.solidTemperature_K);
thermalWork=averageSolidT.*strainEntropy;
dthermalWork=averageSolidT.*(model.thermalLoad'*d.u);
dthermalWork(:,its)=dthermalWork(:,its)+0.5*diag(strainEntropy);
bathHeat=model.bathConductance.*(solidT-bc.solidBathTemperature_K);
solidResidual=model.capacity.*(solidT-s0.solidTemperature_K)+thermalWork+ ...
    dt*(model.thermalStiffness*(solidT-referenceT)+solidHeat+bathHeat);
solidJac=dthermalWork+dt*dsolidHeat;
solidJac(:,its)=solidJac(:,its)+diag(model.capacity)+ ...
    dt*(model.thermalStiffness+diag(model.bathConductance));
affinity=s.phaseAffinity_Jkg+p.phaseViscosity_Jm2skg2*(s.x(:,3)-s0.x(:,3))/dt;
dAffinity=d.affinity+p.phaseViscosity_Jm2skg2*d.ice/dt;
scaleA=p.latentHeat_Jkg/p.meltingTemperature_K;
trial=s.fraction-affinity/scaleA; phase=affinity/scaleA; phaseJac=dAffinity/scaleA;
lower=trial<=0; upper=trial>=1;
phase(lower)=s.fraction(lower); phaseJac(lower,:)=d.fraction(lower,:);
phase(upper)=s.fraction(upper)-1; phaseJac(upper,:)=d.fraction(upper,:);
ms=w*p.rhoWater*p.aperture0_m; es=ms*p.latentHeat_Jkg;
r=[massResidual./ms;phase;energyResidual./es;solidResidual./model.capacity];
J=[massJac./ms;phaseJac;energyJac./es;solidJac./model.capacity];
if nargout<3, return; end
deltaMechanical=0.5*(s.u+s0.u)'*model.stiffness*deltaU;
deltaSolidEnergy=sum(model.capacity.*(solidT-s0.solidTemperature_K))+ ...
    model.referenceTemperature_K*sum(strainEntropy)+deltaMechanical;
deltaFluidEnergy=sum(w.*(s.internalEnergy_Jm2-s0.internalEnergy_Jm2));
netBoundaryEnergy=-dt*sum(model.divergence*flux(:,2))-dt*sum(bathHeat);
entropyChange=sum(w.*(s.entropy_Jm2K-s0.entropy_Jm2K))+ ...
    sum(model.capacity.*log(solidT./s0.solidTemperature_K))+sum(strainEntropy);
entropyIn=dt*boundaryEntropy-dt*sum(bathHeat./bc.solidBathTemperature_K);
mechanicalResidual=model.stiffness*s.u-model.thermalLoad*(solidT-model.referenceTemperature_K)- ...
    model.gapOperator'*(w*p.wetAreaFraction.*s.pressure);
a=struct('state',s,'faceFlux',flux,'localMassDefect_kgm',massResidual, ...
    'localFluidEnergyDefect_Jm',energyResidual,'localSolidEnergyDefect_Jm',solidResidual, ...
    'totalEnergyDefect_Jm',deltaSolidEnergy+deltaFluidEnergy-netBoundaryEnergy, ...
    'entropyProduced_JmK',entropyChange-entropyIn,'faceEntropyRate_WmK',entropyFace, ...
    'pressureWork_Jm',sum(w.*pressureWork),'solidEnergyChange_Jm',deltaSolidEnergy, ...
    'fluidEnergyChange_Jm',deltaFluidEnergy,'boundaryEnergy_Jm',netBoundaryEnergy, ...
    'mechanicalResidual_Nm',mechanicalResidual,'jointHeat_Wm',jointHeat,'solidHeat_Wm',solidHeat, ...
    'phaseAffinity_Jkg',affinity,'scaledResidual',norm(r,inf));
end

function [lstate,rstate,dl,dr,m,dm,k,dk]=faceData(l,r,f,s,d,mob,dmob,kap,dkap,model,bc,distance,nz)
if l>0 && r>0
    lstate=[s.pressure(l),s.x(l,4)]; rstate=[s.pressure(r),s.x(r,4)];
    dl=[d.pressure(l,:);d.temperature(l,:)]; dr=[d.pressure(r,:);d.temperature(r,:)];
    [m,dm]=harmonic(mob(l),mob(r),dmob(l,:),dmob(r,:));
    [k,dk]=harmonic(kap(l),kap(r),dkap(l,:),dkap(r,:));
else
    b=model.boundaryID(f); cellID=max(l,r);
    cellState=[s.pressure(cellID),s.x(cellID,4)];
    reservoir=[bc.reservoirPressure_Pa(b),bc.reservoirTemperature_K(b)];
    cellDerivative=[d.pressure(cellID,:);d.temperature(cellID,:)];
    if l==0
        lstate=reservoir; rstate=cellState; dl=zeros(2,nz); dr=cellDerivative;
    else
        lstate=cellState; rstate=reservoir; dl=cellDerivative; dr=zeros(2,nz);
    end
    m=mob(cellID); dm=dmob(cellID,:); k=kap(cellID); dk=dkap(cellID,:);
    if bc.sealedMass(b) || (bc.outletSeepage(b) && lstate(1)<=rstate(1))
        m=0; dm(:)=0;
    end
    if bc.insulatedHeat(b), k=0; dk(:)=0; end
end
m=m/distance; dm=dm/distance; k=k/distance; dk=dk/distance;
end

function [value,derivative]=harmonic(a,b,da,db)
if a+b==0
    value=0; derivative=zeros(size(da));
else
    value=2*a*b/(a+b); derivative=2*(b^2*da+a^2*db)/(a+b)^2;
end
end
