function [flux,J,a] = premelting_ring_flux(s,model,bc,p)
%PREMELTING_RING_FLUX Conservative liquid and heat ports for a wet ice slit.
% Jacobian J is mass flux wrt liquid pressures, at fixed temperatures.
% Two-film conductance and finite path lengths; no imposed refill source.
q=s.slit; e=q.equilibrium; n=numel(model.weight); nf=size(model.faces,1);
vl=e.liquid.volume(:,1); vlp=e.liquid.volumePressure(:,1); h=e.film_m;
mob=q.chemicalMobility;
dmob=h.^2/(2*p.viscosity_Pas).*e.filmDerivative(:,1)./vl.^2-2*mob.*vlp./vl;
kap=model.conductivityWater*q.phaseVolume_m(:,1)+model.conductivityIce*q.phaseVolume_m(:,2);
flux=zeros(nf,2); J=zeros(nf,n); entropy=zeros(nf,1);
for k=1:nf
    l=model.faces(k,1); r=model.faces(k,2); L=model.faceDistance(k);
    if l>0 && r>0
        left=[e.liquidPressure_Pa(l),e.temperature_K(l)];
        right=[e.liquidPressure_Pa(r),e.temperature_K(r)];
        mobSum=mob(l)+mob(r);
        mf=2*mob(l)*mob(r)/mobSum;
        dmf=zeros(1,n); dmf(l)=2*mob(r)^2/mobSum^2*dmob(l); dmf(r)=2*mob(l)^2/mobSum^2*dmob(r);
        kf=2*kap(l)*kap(r)/(kap(l)+kap(r));
    else
        b=model.boundaryID(k); j=max(l,r);
        cellState=[e.liquidPressure_Pa(j),e.temperature_K(j)];
        reservoir=[bc.reservoirPressure_Pa(b),bc.reservoirTemperature_K(b)];
        if l==0, left=reservoir; right=cellState; else, left=cellState; right=reservoir; end
        mf=mob(j); dmf=zeros(1,n); dmf(j)=dmob(j); kf=kap(j);
        if bc.sealedMass(b) || (bc.outletSeepage(b)&&left(1)<=right(1)), mf=0; dmf(:)=0; end
        if bc.insulatedHeat(b), kf=0; end
    end
    [value,d,aFace]=joint_nonisothermal_face(left,right,mf/L,kf/L,p);
    flux(k,:)=value'; J(k,:)=d(1,5)*dmf/L;
    if l>0, J(k,l)=J(k,l)+d(1,1); end
    if r>0, J(k,r)=J(k,r)+d(1,3); end
    entropy(k)=aFace.entropyProduction_WmK;
end
a=struct('faceEntropyRate_WmK',entropy);
end
