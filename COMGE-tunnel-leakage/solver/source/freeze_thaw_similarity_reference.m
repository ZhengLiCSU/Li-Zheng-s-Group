function reference=freeze_thaw_similarity_reference(surfaceT,farT)
%FREEZE_THAW_SIMILARITY_REFERENCE Independent nonlinear phase-front solution.
% Saturated equal-density soil, constant k=2 W/(m K), no water movement.
% The reference solves C(T)*T_t=k*T_zz in a semi-infinite column using
% eta=z/(2*sqrt(k/Cwarm*t)). Its ODE is T''+2*eta*C(T)/Cwarm*T'=0.
% Capacity is independently differentiated from the published retention
% and enthalpy formulas, without calling the production constitutive law
% or spatial assembly. Both freezing and thawing have a resolved mushy zone.
% This is a numerical similarity reference for that closure, not an exact
% sharp-interface Stefan solution or calibration against soil experiments.
assert(isscalar(surfaceT)&&isscalar(farT)&&isfinite(surfaceT)&&isfinite(farT) ...
    &&surfaceT*farT<0,'PhaseReference:Temperature','Opposite-sign temperatures required.');
p=struct('porosity',.35,'thetaR',.05,'vgAlpha',.02,'vgN',1.67, ...
    'rho',1000,'rhoS',2650,'cs',800,'cw',4180,'ci',2100, ...
    'L',3.34e5,'g',9.81,'Tm',273.15,'conductivity',2);
Cwarm=(1-p.porosity)*p.rhoS*p.cs+p.rho*p.porosity*p.cw;
alpha=p.conductivity/Cwarm;
ode=@(eta,Y)[Y(2);-2*eta*capacity(Y(1),p)/Cwarm*Y(2)];
boundary=@(a,b)[a(1)-surfaceT;b(1)-farT];
guess=@(eta)[surfaceT+(farT-surfaceT)*erf(eta); ...
    (farT-surfaceT)*2/sqrt(pi)*exp(-eta^2)];
initial=bvpinit(linspace(0,8,161),guess);
settings=bvpset('RelTol',1e-8,'AbsTol',1e-10,'NMax',12000);
sol=bvp4c(ode,boundary,initial,settings);
front=fzero(@(eta)deval(sol,eta,1),[0 8]);
% A different collocation order and tighter reference tolerance checks the
% mathematical reference; this does not alter production solver tolerances.
alternative=bvp5c(ode,boundary,bvpinit(sol.x,@(eta)deval(sol,eta)), ...
    bvpset('RelTol',1e-9,'AbsTol',1e-11,'NMax',20000));
checkEta=linspace(0,8,2001);
alternativeTemperatureDifference=max(abs(deval(sol,checkEta,1)-deval(alternative,checkEta,1)));
alternativeFront=fzero(@(eta)deval(alternative,eta,1),[0 8]);
eta=linspace(0,8,20001); [fields,derivative]=deval(sol,eta);
integral= cumtrapz(eta,2*eta.*capacity(fields(1,:),p)/Cwarm);
firstIntegral=fields(2,:).*exp(min(integral,700));
mask=abs(fields(2,:))>1e-9;
firstIntegralRelativeDrift=max(abs(firstIntegral(mask)/firstIntegral(1)-1));
reference=struct('surfaceTemperature_C',surfaceT,'farTemperature_C',farT, ...
    'parameters',p,'warmCapacity_J_m3_K',Cwarm,'alpha_m2_s',alpha, ...
    'solution',sol,'zeroTemperatureEta',front, ...
    'farGradient_dT_dEta',fields(2,end), ...
    'firstIntegralRelativeDrift',firstIntegralRelativeDrift, ...
    'alternativeSolverTemperatureDifference_C',alternativeTemperatureDifference, ...
    'alternativeSolverFrontEtaDifference',abs(front-alternativeFront), ...
    'maxRelativeODEResidual',max(abs(derivative(2,:)+ ...
        2*eta.*capacity(fields(1,:),p)/Cwarm.*fields(2,:)))/max(abs(derivative(2,:))), ...
    'scope','Independent similarity ODE for reversible saturated phase change; not experimental validation.');
end

function C=capacity(T,p)
% Independent direct formula; temperature range is moderate and finite.
sr=p.thetaR/p.porosity;m=1-1/p.vgN;
v=p.vgAlpha*p.L/(p.g*p.Tm)*max(-T,0);
liquid=sr+(1-sr)*(1+v.^p.vgN).^(-m);
liquidT=zeros(size(T));cold=T<0;
liquidT(cold)=(1-sr)*m*p.vgN*p.vgAlpha*p.L/(p.g*p.Tm) ...
    *v(cold).^(p.vgN-1).*(1+v(cold).^p.vgN).^(-m-1);
C=(1-p.porosity)*p.rhoS*p.cs+p.rho*p.porosity*( ...
    p.ci+(p.cw-p.ci)*liquid+(p.L+(p.cw-p.ci)*T).*liquidT);
end
