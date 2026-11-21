function a=soil_poroelastic_fem(xy,tri,T,h,referenceT,settings)
%SOIL_POROELASTIC_FEM P1 plane-strain thermo-poroelastic soil operators.
% A=SOIL_POROELASTIC_FEM(XY,TRI,T,H,TREF,SETTINGS) returns unconstrained
% operators. Displacements are [ux1;uz1;ux2;uz2;...], Z positive down;
% stresses/strains use tension positive and engineering shear gamma_xz.
% SETTINGS fields (all optional except phaseParameters):
% E [Pa], nu, alphaThermal [1/K], biot, rhoBulk [kg/m3], gravity [gx;gz],
% phasePressure ('bishopIce', default, or 'liquid'), phaseParameters (TH p).
% E/nu/alphaThermal/biot/rhoBulk may be scalar or one value per element.
% Body density is prescribed; fbody has no T/h tangent. It is not a fitted
% frost-heave law. Root solver supplies supports, tunnel contact and loads.
% sigma_xz=C_plane*epsilon_xz-(3*Kbulk*alphaThermal*dT+biot*pbar)*[1;1;0].
% The full 3D isotropic thermal eigenstress is used despite epsilon_yy=0.
% Bishop simplification pbar=sl*pw+si*pi (air gauge pressure zero), with
% equal-density Clapeyron candidate pi=pw-rho*L/Tm*min(T,0). This pressure
% participates only where si>0. It is an equilibrium pore-ice candidate,
% not an independently solved ice-lens pressure or lens-initiation model.
% 'liquid' instead takes pbar=pw, useful for a saturated classical Biot
% limit. Neither mode adds a prescribed 9% frost eigenstrain or damage.
% The TH closure is Painter-Karra beta=omega=1, equal rho, fixed solid/pore
% geometry before poroelastic coupling. Mass coupling is assembled by the
% caller from biot*epsilon_v and phase water storage on the SAME SCVs.
% A.V maps u to element div(u); A.nodeProjection is the area/3-volume
% projection of element constants, not an alternative mechanics stencil.
% A.dloadT/H are analytic tangents of fthermal+fpore at fixed referenceT.

if nargin<6||isempty(settings), settings=struct(); end
assert(isnumeric(xy)&&isreal(xy)&&size(xy,2)==2&&all(isfinite(xy(:))), ...
    'soil_poroelastic_fem:InvalidMesh','xy must be finite real ns-by-2 coordinates.');
ns=size(xy,1); ne=size(tri,1);
assert(isnumeric(tri)&&size(tri,2)==3&&isreal(tri)&&all(isfinite(tri(:)))&& ...
    all(tri(:)==fix(tri(:)))&&all(tri(:)>=1&tri(:)<=ns)&&ne>0&&ns>=3, ...
    'soil_poroelastic_fem:InvalidMesh','tri must contain valid three-node triangles.');
T=stateVector(T,ns,'Temperature'); h=stateVector(h,ns,'Head');
referenceT=stateVector(referenceT,ns,'ReferenceTemperature');
assert(isstruct(settings)&&isscalar(settings)&&isfield(settings,'phaseParameters'), ...
    'soil_poroelastic_fem:InvalidSettings','Supply settings.phaseParameters with the TH soil parameters.');
s=struct('E',80e6,'nu',.35,'alphaThermal',1.2e-5,'biot',.65, ...
    'rhoBulk',2000,'gravity',[0;9.81],'phasePressure','bishopIce', ...
    'phaseParameters',settings.phaseParameters);
names=fieldnames(settings);
for j=1:numel(names)
    assert(isfield(s,names{j}),'soil_poroelastic_fem:InvalidSettings', ...
        'Unknown mechanics setting: %s.',names{j});
    s.(names{j})=settings.(names{j});
end
E=elementField(s.E,ne,'E'); nu=elementField(s.nu,ne,'nu');
alpha=elementField(s.alphaThermal,ne,'alphaThermal');
biot=elementField(s.biot,ne,'biot'); density=elementField(s.rhoBulk,ne,'rhoBulk');
assert(all(E>0)&all(nu>-1&nu<.5)&all(biot>=0&biot<=1)&all(density>=0), ...
    'soil_poroelastic_fem:InvalidSettings', ...
    'Require E>0, -1<nu<0.5, 0<=biot<=1, and rhoBulk>=0.');
assert(isnumeric(s.gravity)&&isreal(s.gravity)&&numel(s.gravity)==2&&all(isfinite(s.gravity(:))), ...
    'soil_poroelastic_fem:InvalidSettings','gravity must contain [gx;gz] in m/s2.');
gravity=s.gravity(:); mode=char(s.phasePressure);
assert(any(strcmp(mode,{'liquid','bishopIce'})), ...
    'soil_poroelastic_fem:InvalidSettings','phasePressure must be liquid or bishopIce.');

N=[2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3]; w=ones(3,1)/3;
Tq=T(tri)*N'; hq=h(tri)*N';
[pbar,pT,pH,phase]=phasePressure(Tq,hq,s.phaseParameters,mode);
pMean=pbar*w; pMeanT=pT.*w'; pMeanH=pH.*w';
pMeanT=pMeanT*N; pMeanH=pMeanH*N;
deltaMean=mean(T(tri)-referenceT(tri),2);
nd=2*ns; K=sparse(nd,nd); V=sparse(ne,nd);
A_pressure=sparse(nd,ne); A_temperature=sparse(nd,ne);
fbody=zeros(nd,1); fthermal=fbody; fpore=fbody;
dthermalT=sparse(nd,ns); dporeT=dthermalT; dporeH=dthermalT;
area=zeros(ne,1); grad=zeros(ne,2,3); D=zeros(3,3,ne); B=zeros(3,6,ne);
hydrostatic=[1;1;0]; thermalCoefficient=E.*alpha./(1-2*nu);
for e=1:ne
    ids=tri(e,:); vertices=xy(ids,:);
    cross=(vertices(2,1)-vertices(1,1))*(vertices(3,2)-vertices(1,2))- ...
        (vertices(3,1)-vertices(1,1))*(vertices(2,2)-vertices(1,2));
    area(e)=abs(cross)/2;
    assert(isfinite(area(e))&&area(e)>0, ...
        'soil_poroelastic_fem:InvalidMesh','A triangle has zero area.');
    g=[vertices(2,2)-vertices(3,2),vertices(3,2)-vertices(1,2),vertices(1,2)-vertices(2,2); ...
        vertices(3,1)-vertices(2,1),vertices(1,1)-vertices(3,1),vertices(2,1)-vertices(1,1)]/cross;
    grad(e,:,:)=reshape(g,1,2,3);
    be=zeros(3,6); be(1,1:2:6)=g(1,:); be(2,2:2:6)=g(2,:);
    be(3,1:2:6)=g(2,:); be(3,2:2:6)=g(1,:); B(:,:,e)=be;
    de=E(e)/((1+nu(e))*(1-2*nu(e)))* ...
        [1-nu(e) nu(e) 0;nu(e) 1-nu(e) 0;0 0 (.5-nu(e))];
    D(:,:,e)=de; dof=reshape([2*ids-1;2*ids],[],1);
    K(dof,dof)=K(dof,dof)+area(e)*(be'*de*be);
    V(e,dof)=reshape(g,1,6);
    direction=area(e)*be'*hydrostatic;
    A_pressure(dof,e)=direction*biot(e);
    A_temperature(dof,e)=direction*thermalCoefficient(e);
    fbody(dof)=fbody(dof)+area(e)/3*repmat(density(e)*gravity,3,1);
    fthermal(dof)=fthermal(dof)+direction*(thermalCoefficient(e)*deltaMean(e));
    fpore(dof)=fpore(dof)+direction*(biot(e)*pMean(e));
    dthermalT(dof,ids)=dthermalT(dof,ids)+direction*(thermalCoefficient(e)/3*ones(1,3));
    dporeT(dof,ids)=dporeT(dof,ids)+direction*(biot(e)*pMeanT(e,:));
    dporeH(dof,ids)=dporeH(dof,ids)+direction*(biot(e)*pMeanH(e,:));
end
deposit=sparse(tri(:),repmat((1:ne)',3,1),repmat(area/3,3,1),ns,ne);
volume=full(sum(deposit,2));
assert(all(volume>0),'soil_poroelastic_fem:InvalidMesh','Every soil node must belong to a triangle.');
projection=spdiags(1./volume,0,ns,ns)*deposit;
qrow=repmat((1:ne)',3,1); qcol=tri(:); qvalue=zeros(3*ne,3);
for q=1:3
    values=repmat(N(q,:),ne,1); qvalue(:,q)=values(:);
end
P=sparse([qrow;qrow+ne;qrow+2*ne],repmat(qcol,3,1),qvalue(:),3*ne,ns);
elementAverage=sparse(repmat((1:ne)',3,1),(1:3*ne)',1/3,ne,3*ne);
a=struct('K',K,'fbody',fbody,'fthermal',fthermal,'fpore',fpore, ...
    'dloadT',dthermalT+dporeT,'dloadH',dporeH, ...
    'dthermalT',dthermalT,'dporeT',dporeT,'dporeH',dporeH, ...
    'V',V,'nodeProjection',projection,'nodalVolumetricStrainMatrix',projection*V, ...
    'nodeVolumes',volume,'elementAreas',area,'shapeGradients',grad, ...
    'D',D,'B',B,'thermalEigenstressCoefficient_Pa_K',thermalCoefficient, ...
    'biotElement',biot,'settings',s, ...
    'A_pressure',A_pressure,'A_temperature',A_temperature, ...
    'elementMeanNodal',sparse(repmat((1:ne)',3,1),tri(:),1/3,ne,ns), ...
    'quadrature',struct('shape',N,'weight',w,'P',P,'elementAverage',elementAverage,'T',Tq,'h',hq, ...
        'pw_Pa',phase.pw,'piCandidate_Pa',phase.pi,'pbar_Pa',pbar, ...
        'dpbarT_Pa_K',pT,'dpbarH_Pa_m',pH,'sl',phase.sl,'si',phase.si), ...
    'scope',['Elastic pore-ice pressure coupling; no segregated ice-lens formation, ' ...
        'plastic thaw settlement, cyclic damage, gasket contact or boundary supports.']);
end

function [bar,dT,dH,info]=phasePressure(T,h,p,mode)
phase=frozen_soil_constitutive(T,h,p); rhoG=p.rho*p.g;
pw=rhoG*h; coefficient=p.rho*p.L/p.Tm;
jump=-coefficient*min(T,0); pi=pw+jump;
if strcmp(mode,'liquid')
    bar=pw; dT=zeros(size(T)); dH=rhoG*ones(size(T));
else
    sr=p.thetaR/p.porosity; m=1-1/p.vgN;
    [r,rH]=retention(h,p.vgAlpha,p.vgN,m,sr);
    freezeHead=p.L/(p.g*p.Tm)*min(T,0);
    [F,FH]=retention(freezeHead,p.vgAlpha,p.vgN,m,sr);
    active=F<r;
    slT=zeros(size(T)); slT(active)=FH(active)*(p.L/(p.g*p.Tm));
    slH=rH; slH(active)=0;
    siT=zeros(size(T)); siT(active)=-slT(active)./r(active);
    siH=zeros(size(T)); siH(active)=phase.sl(active).*rH(active)./r(active).^2;
    bar=(phase.sl+phase.si).*pw+phase.si.*jump;
    dT=(slT+siT).*pw+siT.*jump-phase.si*coefficient.*double(T<0);
    dH=(slH+siH).*pw+(phase.sl+phase.si)*rhoG+siH.*jump;
end
info=struct('pw',pw,'pi',pi,'sl',phase.sl,'si',phase.si);
end

function [s,ds]=retention(h,alpha,n,m,sr)
s=ones(size(h)); ds=zeros(size(h)); negative=h<0;
u=-h(negative); lx=n*(log(alpha)+log(u)); l1=max(lx,0)+log1p(exp(-abs(lx)));
lq=-m*l1; s(negative)=sr+(1-sr)*exp(lq);
ds(negative)=exp(log((1-sr)*m*n)+lq+lx-l1-log(u));
end

function v=stateVector(v,n,name)
assert(isnumeric(v)&&isreal(v)&&all(isfinite(v(:)))&& ...
    (isscalar(v)||isvector(v)&&numel(v)==n), ...
    'soil_poroelastic_fem:InvalidState','%s must be scalar or contain one finite value per node.',name);
if isscalar(v), v=repmat(double(v),n,1); else, v=double(v(:)); end
end

function v=elementField(v,n,name)
assert(isnumeric(v)&&isreal(v)&&all(isfinite(v(:)))&& ...
    (isscalar(v)||isvector(v)&&numel(v)==n), ...
    'soil_poroelastic_fem:InvalidSettings','%s must be scalar or contain one finite value per element.',name);
if isscalar(v), v=repmat(double(v),n,1); else, v=double(v(:)); end
end
