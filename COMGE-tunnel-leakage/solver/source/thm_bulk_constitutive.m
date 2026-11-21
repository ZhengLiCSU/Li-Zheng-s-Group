function a=thm_bulk_constitutive(T,h,epsilon,p,existingFrozenState,wantDerivatives)
%THM_BULK_CONSTITUTIVE Current equal-density THM bulk law and exact partials.
% A=thm_bulk_constitutive(T,H,EPSILON,P) applies the v16 reference-volume
% coupling to frozen_soil_constitutive at arbitrary SCV/nodal points.
% Optional fifth input is the UNMODIFIED frozen_soil_constitutive(T,H,P)
% output at the SAME points and parameters, avoiding a second evaluation.
% The caller must not reuse an already-deformed THM state as this input.
% Optional sixth logical input defaults to true; false computes identical
% properties without tangents, accepting a properties-only fifth input.
%
% EPSILON=div(u-uInstalled) is a scalar or an array matching T/H. Retention
% and phase fractions keep the original phi0-based Painter-Karra shape.
% Vp=phi0+alphaB*epsilon, J=1+epsilon, phiCurrent=Vp/J. Storage uses Vp per
% REFERENCE bulk volume; conductivity uses the CURRENT mixture fraction.
% The fixed reference-domain transport geometry is unchanged by this helper.
%
% Only mass/energy/k/K and their T/H partials are overwritten. Original
% thetaW/I/A/total retain the frozen-kernel phi0 convention; currentThetaW/I/A
% separately report current-volume fractions. Phase derivatives slT/slH/
% siT/siH are unchanged. dMassEps,dEnergyEps,dkEps,dKEps are partials with
% respect to volumetric strain, at fixed T/H and the installed reference.
%
% This is precisely the previous v16 constitutive block, including optional
% Kozeny porosity feedback, equal phase densities and arbitrary enthalpy
% datum. No finite deformation, density-expansion eigenstrain or damage law
% is added. The caller still enforces full small-strain/gradient bounds.

assert(nargin>=4&&nargin<=6,'THMBulk:InputCount','Supply T,H,epsilon,P and optional frozen state/derivative switch.');
if nargin<6, wantDerivatives=true; end
assert(islogical(wantDerivatives)&&isscalar(wantDerivatives), ...
    'THMBulk:DerivativeSwitch','wantDerivatives must be a scalar logical.');
if nargin<5||isempty(existingFrozenState)
    a=frozen_soil_constitutive(T,h,p,wantDerivatives);
else
    a=existingFrozenState;
    required={'mass','energy','k','K','sl','si','sg'};
    if wantDerivatives
        required=[required,{'slT','slH','siT','siH', ...
            'dMassT','dMassH','dEnergyT','dEnergyH','dkT','dkH','dKT','dKH'}];
    end
    assert(isstruct(a)&&isscalar(a)&&all(isfield(a,required)), ...
        'THMBulk:FrozenState','Supply an unmodified frozen-kernel state with requested properties/partials.');
    if wantDerivatives
        for j=8:numel(required)
            assert(isequal(size(a.(required{j})),size(a.mass)), ...
                'THMBulk:FrozenState','The requested frozen-kernel partials must match the property array.');
        end
    end
end
assert(isstruct(p)&&isscalar(p)&&isfield(p,'mechanics') ...
    &&isstruct(p.mechanics)&&isscalar(p.mechanics)&&isfield(p.mechanics,'biot'), ...
    'THMBulk:Parameters','P.mechanics.biot is required.');
alpha=p.mechanics.biot;
assert(isnumeric(alpha)&&isreal(alpha)&&isscalar(alpha)&&isfinite(alpha)&&alpha>=0&&alpha<=1, ...
    'THMBulk:Parameters','Biot coefficient must be a finite scalar in [0,1].');
shape=size(a.mass);
T=pointValue(T,shape,'Temperature'); h=pointValue(h,shape,'Head');
epsilon=pointValue(epsilon,shape,'Strain');
phi0=p.porosity; Vp=phi0+alpha*epsilon; J=1+epsilon;
assert(all(J(:)>0)&&all(Vp(:)>0&Vp(:)<J(:)), ...
    'THMBulk:PoreVolumeDomain','Require positive Jacobian and 0<reference pore volume<J.');
phi=Vp./J; deltaVp=Vp-phi0;
reference=0; if isfield(p,'enthalpyReference_Jkg'), reference=p.enthalpyReference_Jkg; end
liquid=p.L+reference+p.cw*T; ice=reference+p.ci*T;
total=a.sl+a.si;
phaseEnthalpy=a.sl.*liquid+a.si.*ice;
a.mass=a.mass+deltaVp.*total;
a.energy=a.energy+p.rho*deltaVp.*phaseEnthalpy;
phaseConductivity=p.kw*a.sl+p.ki*a.si+p.ka*a.sg;
a.k=(1-phi)*p.ks+phi.*phaseConductivity;
feedback=false;
if isfield(p,'thm')&&isfield(p.thm,'permeabilityPorosityFeedback')
    feedback=p.thm.permeabilityPorosityFeedback;
    assert((isnumeric(feedback)||islogical(feedback))&&isreal(feedback) ...
        &&isscalar(feedback)&&isfinite(feedback)&&ismember(feedback,[0 1]), ...
        'THMBulk:Parameters','permeabilityPorosityFeedback must be scalar logical 0/1.');
end
if feedback
    multiplier=(phi/phi0).^3.*((1-phi0)./(1-phi)).^2;
    a.K=a.K.*multiplier;
end
a.referencePoreVolume=Vp; a.currentPorosity=phi;
a.referenceJacobian=J; a.volumetricStrain=epsilon;
a.currentThetaW=phi.*a.sl; a.currentThetaI=phi.*a.si; a.currentThetaA=phi.*a.sg;
if ~wantDerivatives
    derivatives={'slT','slH','siT','siH','dMassT','dMassH', ...
        'dEnergyT','dEnergyH','dkT','dkH','dKT','dKH', ...
        'dMassEps','dEnergyEps','dkEps','dKEps', ...
        'dReferencePoreVolumeEps','dCurrentPorosityEps'};
    for j=1:numel(derivatives), a.(derivatives{j})=[]; end
    return
end
dPhi=(alpha-phi0)./J.^2;
totalT=a.slT+a.siT; totalH=a.slH+a.siH;
phaseEnthalpyT=a.slT.*liquid+a.siT.*ice+p.cw*a.sl+p.ci*a.si;
phaseEnthalpyH=a.slH.*liquid+a.siH.*ice;
a.dMassT=a.dMassT+deltaVp.*totalT;
a.dMassH=a.dMassH+deltaVp.*totalH;
a.dMassEps=alpha*total;
a.dEnergyT=a.dEnergyT+p.rho*deltaVp.*phaseEnthalpyT;
a.dEnergyH=a.dEnergyH+p.rho*deltaVp.*phaseEnthalpyH;
a.dEnergyEps=p.rho*alpha*phaseEnthalpy;
a.dkT=phi.*((p.kw-p.ka)*a.slT+(p.ki-p.ka)*a.siT);
a.dkH=phi.*((p.kw-p.ka)*a.slH+(p.ki-p.ka)*a.siH);
a.dkEps=(phaseConductivity-p.ks).*dPhi;
a.dKEps=zeros(shape);
if feedback
    a.dKT=a.dKT.*multiplier; a.dKH=a.dKH.*multiplier;
    a.dKEps=a.K.*(3./phi+2./(1-phi)).*dPhi;
end
a.dReferencePoreVolumeEps=alpha*ones(shape); a.dCurrentPorosityEps=dPhi;
end

function value=pointValue(value,shape,name)
assert(isnumeric(value)&&isreal(value)&&all(isfinite(value(:))) ...
    &&(isscalar(value)||isequal(size(value),shape)), ...
    ['THMBulk:' name],'%s must be finite scalar or match the constitutive point array.',name);
value=double(value); if isscalar(value), value=repmat(value,shape); end
end
