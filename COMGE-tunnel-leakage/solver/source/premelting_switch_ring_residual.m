function [r,J,a] = premelting_switch_ring_residual(z,old,dt,model,bc,p,f,history,loadingAtHistory,oldState)
%PREMELTING_SWITCH_RING_RESIDUAL Conserved three-domain premelting THM branch.
% Independent joint T, liquid chemical pressure, three domain areas and
% wall pressure; host temperatures are unknowns. Midpoint pressure work and reciprocal thermoelastic heating
% make TOTAL internal energy conservative; no reconstructed thermostat heat.
% Nonlinear mechanical equilibrium and history-dependent tangent are condensed exactly.
if nargin<9, loadingAtHistory=[]; end
needJacobian=nargout>1;
if needJacobian
    [s,d]=premelting_switch_ring_state(z,model,p,f,history,loadingAtHistory);
else
    s=premelting_switch_ring_state(z,model,p,f,history,loadingAtHistory);
end
% The step supplies its immutable old-state cache; uncached calls are valid.
if nargin<10, s0=premelting_switch_ring_state(old,model,p,f,history); else, s0=oldState; end
n=numel(model.weight); ns=numel(model.capacity); nz=numel(z); w=model.weight;
its=6*n+(1:ns);
vw=1/p.rhoWater; cw=vw/p.bulkWater_Pa;
if ~needJacobian
    nz=0; empty=zeros(n,0);
    d=struct('chemicalMobility',empty,'liquidVolume',empty,'iceVolume',empty, ...
        'vaporVolume',empty,'liquidPressure',empty,'temperature',empty);
end
mob=s.chemicalMobility; dmob=d.chemicalMobility;
kappa=model.conductivityWater*s.liquidVolume_m+model.conductivityIce*s.iceVolume_m+p.vaporConductivity_WmK*s.vaporVolume_m;
dkappa=model.conductivityWater*d.liquidVolume+model.conductivityIce*d.iceVolume+p.vaporConductivity_WmK*d.vaporVolume;
nf=size(model.faces,1); flux=zeros(nf,2); dmass=zeros(nf,nz); denergy=dmass;
entropyFace=zeros(nf,1); boundaryEntropy=0;
for f=1:nf
    l=model.faces(f,1); q=model.faces(f,2); distance=model.faceDistance(f);
    [left,right,dl,dr,m,dm,k,dk]=faceData(l,q,f,s,d,mob,dmob,kappa,dkappa,model,bc,distance,nz);
    if needJacobian
        [value,derivative,e]=joint_nonisothermal_face(left,right,m,k,p);
        chain=[dl;dr;dm;dk]; gradient=derivative*chain;
        dmass(f,:)=gradient(1,:); denergy(f,:)=gradient(2,:);
        entropyFace(f)=e.entropyProduction_WmK;
    else
        value=joint_nonisothermal_face(left,right,m,k,p);
    end
    flux(f,:)=value';
    if nargout>2 && (l==0 || q==0)
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
if needJacobian
    djointHeat=model.exchangeRowSum.*d.temperature; djointHeat(:,its)=djointHeat(:,its)-H;
    dsolidHeat=-H'*d.temperature; dsolidHeat(:,its)=dsolidHeat(:,its)+diag(model.exchangeColumnSum);
end
deltaGap=s.x(:,1)-s0.x(:,1); pAverage=0.5*(s.pressure+s0.pressure);
pressureWork=p.wetAreaFraction*pAverage.*deltaGap;
massResidual=w.*(s.x(:,2)-s0.x(:,2))+dt*model.divergence*flux(:,1);
energyResidual=w.*(s.internalEnergy_Jm2-s0.internalEnergy_Jm2+pressureWork)+ ...
    dt*(model.divergence*flux(:,2)+jointHeat);
if needJacobian
    dwork=p.wetAreaFraction*(0.5*deltaGap.*d.pressure+pAverage.*d.gap);
    massJac=w.*d.mass+dt*model.divergence*dmass;
    energyJac=w.*(d.energy+dwork)+dt*(model.divergence*denergy+djointHeat);
end
% Reciprocal term T d(thermalLoadMap' U) is REQUIRED by the solid free energy.
deltaU=s.u-s0.u; strainEntropy=model.thermalLoad'*deltaU;
averageSolidT=0.5*(solidT+s0.solidTemperature_K);
thermalWork=averageSolidT.*strainEntropy;
bathHeat=model.bathConductance.*(solidT-bc.solidBathTemperature_K);
solidResidual=model.capacity.*(solidT-s0.solidTemperature_K)+thermalWork+ ...
    dt*(model.thermalStiffness*(solidT-referenceT)+solidHeat+bathHeat);
scaleA=p.rhoWater*p.aperture0_m*p.latentHeat_Jkg/p.meltingTemperature_K;
wallResidual=(s.pressure-s.wallPressure_Pa)/1e6;
if needJacobian
    dthermalWork=averageSolidT.*(model.thermalLoad'*d.u);
    dthermalWork(:,its)=dthermalWork(:,its)+0.5*diag(strainEntropy);
    solidJac=dthermalWork+dt*dsolidHeat;
    solidJac(:,its)=solidJac(:,its)+diag(model.capacity)+ ...
        dt*(model.thermalStiffness+diag(model.bathConductance));
    [phase,pt,po]=phase_simplex_residual(s.areaFraction,s.grandPotential_Jm2,scaleA);
    phaseJac=zeros(3*n,nz);
    for cellID=1:n
        localJ=reshape(pt(cellID,:,:),3,3)*reshape(d.areas(cellID,:,:),3,nz)+ ...
            reshape(po(cellID,:,:),3,3)*reshape(d.grandPotential(cellID,:,:),3,nz);
        phaseJac(cellID+n*(0:2),:)=localJ;
    end
    wallJac=(d.pressure-d.wallPressure)/1e6;
else
    phase=phase_simplex_residual(s.areaFraction,s.grandPotential_Jm2,scaleA);
end
ms=w.*p.rhoWater.*p.aperture0_m; es=ms*p.latentHeat_Jkg;
r=[massResidual./ms;phase(:);energyResidual./es;wallResidual;solidResidual./model.capacity];
if ~needJacobian, return; end
J=[massJac./ms;phaseJac;energyJac./es;wallJac;solidJac./model.capacity];
if nargout<3, return; end
deltaMechanical=0.5*(s.u+s0.u)'*model.stiffness*deltaU;
deltaElastic=s.cohesive.elasticEnergy_Jm2-s0.cohesive.elasticEnergy_Jm2;
deltaSurface=s.cohesive.surfaceEnergy_Jm2-s0.cohesive.surfaceEnergy_Jm2;
cohesiveWork=0.5*(s.cohesive.traction_Pa+s0.cohesive.traction_Pa).*deltaGap;
% Fracture surface energy is explicitly retained in total internal energy.
% Midpoint-work quadrature is exact on each linear loading/unloading piece.
% Crossing a kink produces a REPORTED defect and triggers step refinement;
% it is never injected into temperature as artificial heat.
quadratureDefect=deltaElastic+deltaSurface-cohesiveWork;
deltaSolidEnergy=sum(w.*(deltaElastic+deltaSurface))+sum(model.capacity.*(solidT-s0.solidTemperature_K))+ ...
    model.referenceTemperature_K*sum(strainEntropy)+deltaMechanical;
deltaFluidEnergy=sum(w.*(s.internalEnergy_Jm2-s0.internalEnergy_Jm2));
netBoundaryEnergy=-dt*sum(model.divergence*flux(:,2))-dt*sum(bathHeat);
entropyChange=sum(w.*(s.entropy_Jm2K-s0.entropy_Jm2K))+ ...
    sum(model.capacity.*log(solidT./s0.solidTemperature_K))+sum(strainEntropy);
entropyIn=dt*boundaryEntropy-dt*sum(bathHeat./bc.solidBathTemperature_K);
mechanicalResidual=model.stiffness*s.u+model.gapOperator'*(w.*s.cohesive.traction_Pa)-model.thermalLoad*(solidT-model.referenceTemperature_K)- ...
    model.gapOperator'*(w*p.wetAreaFraction.*s.pressure);
a=struct('state',s,'fractureSurfaceEnergyChange_Jm',sum(w.*deltaSurface), ...
    'cohesiveQuadratureDefect_Jm',sum(w.*quadratureDefect), ...
    'localCohesiveQuadratureDefect_Jm2',quadratureDefect, ...
    'surfaceEnergy_Jm',sum(w.*s.cohesive.surfaceEnergy_Jm2),'faceFlux',flux,'localMassDefect_kgm',massResidual, ...
    'localFluidEnergyDefect_Jm',energyResidual,'localSolidEnergyDefect_Jm',solidResidual, ...
    'totalEnergyDefect_Jm',deltaSolidEnergy+deltaFluidEnergy-netBoundaryEnergy, ...
    'entropyProduced_JmK',entropyChange-entropyIn,'faceEntropyRate_WmK',entropyFace, ...
    'pressureWork_Jm',sum(w.*pressureWork),'solidEnergyChange_Jm',deltaSolidEnergy, ...
    'fluidEnergyChange_Jm',deltaFluidEnergy,'boundaryEnergy_Jm',netBoundaryEnergy, ...
    'mechanicalResidual_Nm',mechanicalResidual,'jointHeat_Wm',jointHeat,'solidHeat_Wm',solidHeat, ...
    'phaseResidual',phase,'wallResidual_Pa',1e6*wallResidual,'scaledResidual',norm(r,inf));
end

function [lstate,rstate,dl,dr,m,dm,k,dk]=faceData(l,r,f,s,d,mob,dmob,kap,dkap,model,bc,distance,nz)
if l>0 && r>0
    lstate=[s.liquidPressure(l),s.x(l,4)]; rstate=[s.liquidPressure(r),s.x(r,4)];
    dl=[d.liquidPressure(l,:);d.temperature(l,:)]; dr=[d.liquidPressure(r,:);d.temperature(r,:)];
    [m,dm]=harmonic(mob(l),mob(r),dmob(l,:),dmob(r,:));
    [k,dk]=harmonic(kap(l),kap(r),dkap(l,:),dkap(r,:));
else
    b=model.boundaryID(f); cellID=max(l,r);
    cellState=[s.liquidPressure(cellID),s.x(cellID,4)];
    reservoir=[bc.reservoirPressure_Pa(b),bc.reservoirTemperature_K(b)];
    cellDerivative=[d.liquidPressure(cellID,:);d.temperature(cellID,:)];
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
