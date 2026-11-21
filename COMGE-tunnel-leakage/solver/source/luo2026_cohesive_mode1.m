function [traction,tangent,memory]=luo2026_cohesive_mode1(gap,fatigue,old,material)
%LUO2026_COHESIVE_MODE1 Mode-I restriction of Eqs.18,20-25 (beta=1).
% Constitutive reference only; NOT a gasket law in the production tunnel.
% Compressive penalty stays intact. Positive traction is tensile.
arguments
    gap (1,1) double {mustBeReal,mustBeFinite}
    fatigue (1,1) double {mustBeReal,mustBeFinite}
    old (1,1) struct
    material (1,1) struct
end
assert(all(isfield(material,{'strength_Pa','crackingGap_m','debondingGap_m'})), ...
    'Luo2026:Material','Supply strength and cracking/debonding separation.');
s=material.strength_Pa;r0=material.crackingGap_m;rc=material.debondingGap_m;
assert(s>0&&isfinite(s)&&r0>0&&rc>r0&&isfinite(rc)&&fatigue>=0&&fatigue<=1, ...
    'Luo2026:Material','Invalid material or fatigue damage.');
if isempty(fieldnames(old)),old.omegaMax=0;end
assert(isfield(old,'omegaMax')&&isscalar(old.omegaMax)&&isfinite(old.omegaMax) ...
    &&old.omegaMax>=0,'Luo2026:Memory','Invalid irreversible separation memory.');
omega=max(max(gap,0)/r0-1,0);bar=max(old.omegaMax,omega);eta=1-r0/rc;
D=min(bar/(1+bar)/eta,1);K=s/r0;
memory=struct('omegaMax',bar,'separationDamage',D,'fatigueDamage',fatigue);
if gap<0
    traction=K*gap;tangent=K;
else
    traction=(1-D)*(1-fatigue)*K*gap;
    tangent=(1-D)*(1-fatigue)*K;
    if omega>old.omegaMax&&gap>r0&&gap<rc
        tangent=-(1-fatigue)*s/(rc-r0);
    end
end
end
