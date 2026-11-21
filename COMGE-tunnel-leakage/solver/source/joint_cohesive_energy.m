function s = joint_cohesive_energy(gap,kappaOld,c)
%JOINT_COHESIVE_ENERGY Irreversible bilinear cohesive law + intact contact.
% Per reference interface area. History is returned, never committed here.
% Psi=.5*(1-d)*Kt*<g>+^2+.5*Kc*<g>-^2.
% Gamma=Gc*(kappa-g0)/(gf-g0) in softening, bounded by Gc.
% Under active damage Gamma'(d)=.5*Kt*kappa^2: energetic consistency.
arguments
    gap (:,1) double {mustBeFinite,mustBeReal}
    kappaOld (:,1) double {mustBeFinite,mustBeNonnegative}
    c (1,1) struct = joint_cohesive_parameters()
end
assert(numel(gap)==numel(kappaOld),'Cohesive:Shape','History must match the interface.');
g0=c.onsetGap_m; gf=c.failureGap_m; kt=c.tensileStiffness_Pam; kc=c.contactStiffness_Pam;
kappa=max(kappaOld,max(gap,0)); d=zeros(size(gap)); dd=zeros(size(gap));
secant=kt*ones(size(gap)); active=false(size(gap));
surface=zeros(size(gap)); dSurface=zeros(size(gap));
if c.enabled
    soft=kappa>g0 & kappa<gf; broken=kappa>=gf;
    secant(soft)=c.strength_Pa*(gf-kappa(soft))./((gf-g0)*kappa(soft));
    secant(broken)=0; d=1-secant/kt;
    active=soft & gap>kappaOld;
    dd(active)=gf*g0./((gf-g0)*kappa(active).^2);
    surface(soft)=c.fractureEnergy_Jm2*(kappa(soft)-g0)/(gf-g0);
    surface(broken)=c.fractureEnergy_Jm2;
    dSurface(active)=c.fractureEnergy_Jm2/(gf-g0);
end
positive=max(gap,0); negative=min(gap,0);
traction=secant.*positive+kc*negative;
tangent=secant.*(gap>=0)+kc*(gap<0);
tangent(active)=-c.strength_Pa/(gf-g0);
elastic=0.5*secant.*positive.^2+0.5*kc*negative.^2;
s=struct('traction_Pa',traction,'tangent_Pam',tangent,'kappa_m',kappa, ...
    'damage',d,'damageDerivative',dd,'elasticEnergy_Jm2',elastic, ...
    'surfaceEnergy_Jm2',surface,'surfaceDerivative',dSurface, ...
    'elasticDerivative',traction-0.5*kt*positive.^2.*dd, ...
    'tensionReleaseRate_Jm2',0.5*kt*positive.^2);
end
