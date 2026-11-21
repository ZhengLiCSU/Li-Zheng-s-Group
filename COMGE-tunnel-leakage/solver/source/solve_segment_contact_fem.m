function out=solve_segment_contact_fem(M,Tglobal,settings)
%SOLVE_SEGMENT_CONTACT_FEM Standalone eight-block plane-strain Q4 mechanics.
% OUT=SOLVE_SEGMENT_CONTACT_FEM(M,TGLOBAL,SETTINGS) uses a mesh from
% BUILD_SEGMENT_MECHANICAL_MESH. TGLOBAL is the original TH nodal temperature
% in C; SETTINGS.referenceTemperatureGlobal_C is its scalar or nodal reference.
% This is ring-only mechanics with prescribed elastic outer support. It is
% not a soil equilibrium solve, a geostatic/construction analysis, or THM.
% No gasket sealing law, seepage feedback, frost heave or damage is added.
%
% Unilateral normal contact is frictionless and compression-only. One
% discrete bolt per joint/per metre of tunnel axis carries tension only;
% its shear connector is optional. Initial clearances and prestress are
% explicit inputs. Joint gap, slip and relative face rotation come from U.
% Penalty contact permits small overclosure; it does not enforce exact gap=0.
%
% Foundation stiffness [Pa/m] and displacement are prescribed inputs, not
% calibrated soil pressure. Displacement may be supplied at all mechanical
% nodes in foundationDisplacement_m, or at original global SOIL nodes in
% foundationDisplacementGlobal_m via the preserved soil-interface map.
% Fixed DOFs, if needed, are explicit. Uncontrolled rigid modes are rejected.
%
% Set solveEquilibrium=false to evaluate the public residual/tangent at
% initialDisplacement_m, for constitutive and independent directional tests.
% Cached assembly: C=solve_segment_contact_fem('assemble',M,settings).
% Cached evaluation: E=solve_segment_contact_fem('evaluate',C,state), with
% state.temperatureGlobal_C, state.ringDisplacement_m and optional
% state.soilDisplacementGlobal_m (original-global N-by-2). This evaluation
% performs no geometric assembly and is intended for a coupled solver.
% Supplied soil displacement overrides prescribed foundation data; when
% omitted, the foundation target supplied at assemble time is preserved.
% Set state.computeStress=true for Gauss stress/strain output from cached
% quadrature; leave it false during nonlinear residual evaluations.
% All numerical/material defaults below are illustrative and uncalibrated.
if nargin<3, settings=struct(); end
if ischar(M)||isstring(M)
 command=char(M);
 if strcmp(command,'assemble')
  mesh=Tglobal; options=settings;
  options.solveEquilibrium=false; options.initialDisplacement_m=[];
  response=solve_segment_contact_fem(mesh,20*ones(size(mesh.thermalInterpolation,2),1),options);
  out=response.assembly; return
 elseif strcmp(command,'evaluate')
  out=evaluateCached(Tglobal,settings); return
 else
  error('SegmentContactFEM:Command','Use assemble or evaluate for the cached interface.');
 end
end
assert(isstruct(settings)&&isscalar(settings),'SegmentContactFEM:Settings', ...
 'Settings must be a scalar struct.');
p=defaults(); names=fieldnames(settings);
for k=1:numel(names)
 assert(isfield(p,names{k}),'SegmentContactFEM:Settings','Unknown setting: %s.',names{k});
 p.(names{k})=settings.(names{k});
end
required={'p','q','nodeSegmentID','jointNodePairs','jointInterface','thermalInterpolation','center'};
assert(isstruct(M)&&all(isfield(M,required)),'SegmentContactFEM:Mesh', ...
 'A separate eight-block mechanical mesh is required.');
N=size(M.p,1); nd=2*N; nj=size(M.jointNodePairs,3);
assert(nj==8&&numel(unique(M.nodeSegmentID))==8,'SegmentContactFEM:Mesh', ...
 'This solver requires the eight independent ring blocks.');
assert(isnumeric(Tglobal)&&isreal(Tglobal)&&isvector(Tglobal)&& ...
 numel(Tglobal)==size(M.thermalInterpolation,2)&&all(isfinite(Tglobal(:))), ...
 'SegmentContactFEM:Temperature','Supply a finite original global temperature vector.');
T0=p.referenceTemperatureGlobal_C;
if isscalar(T0), T0=repmat(T0,numel(Tglobal),1); end
assert(isvector(T0)&&numel(T0)==numel(Tglobal)&&all(isfinite(T0(:))), ...
 'SegmentContactFEM:Temperature','Reference temperature must be scalar or original global nodal values.');
deltaT=M.thermalInterpolation*(Tglobal(:)-T0(:));
validateSettings(p,nd,nj);
kn=perJoint(p.contactNormalStiffness_Pa_per_m,nj);
clearance=perJoint(p.initialJointGap_m,nj);
kb=perJoint(p.boltNormalStiffness_N_per_m2,nj);
preload=perJoint(p.boltPreload_N_per_m,nj);
ks=perJoint(p.boltShearStiffness_N_per_m2,nj);
position=perJoint(p.boltRadialFraction,nj);
assert(all(kb>0|preload==0),'SegmentContactFEM:Settings', ...
 'A nonzero bolt preload needs a nonzero bolt stiffness.');
[Kelastic,fthermal,quadrature,thermalMap]=elasticAssembly(M,deltaT,p);
[Kfoundation,foundationLoad,foundationConstant,coupling]=foundationAssembly(M,p);
[An,At,weight,jointID]=jointGaussOperators(M);
g0=clearance(jointID); contactK=kn(jointID).*weight;
[Bn,Bt]=boltOperators(M,position);
fixed=p.fixedDofs(:); fixedValues=p.fixedValues_m(:);
if isempty(fixedValues), fixedValues=zeros(numel(fixed),1); end
free=setdiff((1:nd)',fixed);
if p.solveEquilibrium&&isempty(fixed)&&nnz(Kfoundation)==0
 error('SegmentContactFEM:RigidModes', ...
  'An unsupported ring has free global rigid motions; supply physical support or fixed DOFs.');
end
u=displacementVector(p.initialDisplacement_m,N); u(fixed)=fixedValues;
base=struct('K',Kelastic+Kfoundation,'f',fthermal+foundationLoad, ...
 'foundationConstant',foundationConstant,'An',An,'At',At,'Bn',Bn,'Bt',Bt, ...
 'g0',g0,'contactK',contactK,'contactMaterialK',kn(jointID), ...
 'jointID',jointID,'weight',weight,'kb',kb,'preload',preload,'ks',ks);
[r,J,energy]=evaluate(u,base);
forceScale=max([1,norm(fthermal(free),inf),norm(foundationLoad(free),inf),norm(r(free),inf)]);
limit=p.absoluteForceTolerance_N_per_m+p.relativeTolerance*forceScale;
iterations=0; backtracks=0; converged=norm(r(free),inf)<=limit;
if p.solveEquilibrium
 for iteration=1:p.maximumNewtonIterations
  if converged, break; end
  iterations=iteration;
  [~,flag]=chol(J(free,free));
  assert(flag==0,'SegmentContactFEM:RigidModes', ...
   ['Mechanical tangent has an uncontrolled rigid mode or is not positive definite. ' ...
    'Supply physically justified foundation support or explicit fixed DOFs.']);
  du=zeros(nd,1); du(free)=-J(free,free)\r(free);
  directional=r'*du;
  assert(isfinite(directional)&&directional<0,'SegmentContactFEM:Newton', ...
   'Newton direction is not a finite descent direction.');
  accepted=false; alpha=1;
  for trial=0:p.maximumLineSearchReductions
   candidate=u+alpha*du;
   [rt,Jt,et]=evaluate(candidate,base);
   roundoff=64*eps*max([1,abs(energy),abs(et)]);
   if et<=energy+p.armijo*alpha*directional+roundoff
    accepted=true; break
   end
   alpha=alpha/2;
  end
  assert(accepted,'SegmentContactFEM:LineSearch','Contact Newton line search failed.');
  u=candidate; r=rt; J=Jt; energy=et; backtracks=backtracks+trial;
  converged=norm(r(free),inf)<=limit;
 end
 assert(converged,'SegmentContactFEM:Nonconvergence', ...
  'Contact equilibrium did not converge; residual %.6g, limit %.6g N/m.', ...
  norm(r(free),inf),limit);
 % Detect an unanchored equilibrium even when its zero residual is trivial.
 if ~isempty(free)
  [~,flag]=chol(J(free,free));
  assert(flag==0,'SegmentContactFEM:RigidModes', ...
   'A zero residual does not remove uncontrolled mechanical rigid modes.');
 end
end
U=reshape(u,2,[])';
[gap,slip,rotation]=jointNodalState(M,U,clearance);
boltForce=max(preload+kb.*(Bn*u),0); boltShear=ks.*(Bt*u);
gaussGap=g0+An*u; gaussContact=-base.contactMaterialK.*min(gaussGap,0);
nodalPressure=-min(gap,0).*kn';
[stress,strain]=stressOutput(M,U,deltaT,p,quadrature);
reaction=zeros(nd,1); reaction(fixed)=r(fixed);
status='evaluated_supplied_displacement';
if p.solveEquilibrium, status='contact_equilibrium_converged'; end
out=struct('status',status,'U_m',U,'uVector_m',u,'residual_N_per_m',r, ...
 'tangent_N_per_m2',J,'elasticStiffness_N_per_m2',Kelastic, ...
 'foundationStiffness_N_per_m2',Kfoundation,'thermalLoad_N_per_m',fthermal, ...
 'jointGap_m',gap,'jointSlip_m',slip,'relativeFaceRotation_rad',rotation, ...
 'jointContactPressure_Pa',nodalPressure,'jointGapGauss_m',gaussGap, ...
 'jointContactPressureGauss_Pa',gaussContact,'jointGaussID',jointID, ...
 'boltForce_N_per_m',boltForce,'boltShearForce_N_per_m',boltShear, ...
 'stressGauss_Pa',stress,'strainGauss',strain,'stressComponents',{{'xx','yy','xy','zz'}}, ...
 'stressCoordinateSystem', ...
 'x,y are original S.pAll coordinates (y is positive-down depth); z is the out-of-plane tunnel axis.', ...
 'deltaTemperatureMechanical_C',deltaT,'fixedReaction_N_per_m',reaction, ...
 'forceScale_N_per_m',forceScale,'forceTolerance_N_per_m',limit, ...
 'maximumFreeResidual_N_per_m',norm(r(free),inf), ...
 'normalizedFreeResidual',norm(r(free),inf)/forceScale, ...
 'potentialEnergy_N',energy,'newtonIterations',iterations,'lineSearchReductions',backtracks, ...
 'activeContactGaussCount',nnz(gaussGap<0),'activeTensionBoltCount',nnz(boltForce>0), ...
 'parameters',p,'parameterStatus','Illustrative and uncalibrated; foundation is a prescribed support.', ...
 'scope','Standalone ring contact/bolt mechanics; no soil solve, gasket law, damage or hydraulic feedback.');
soilResidual=coupling.KringSoil'*u+coupling.soilTargetLoad;
out.soilInterfaceResidual_N_per_m=soilResidual;
out.soilReactionForce_global_N_per_m=reshape(-soilResidual,2,[])';
out.interfaceCoupling=coupling;
[nodalNormal,nodalTangent]=nodalJointOperators(M);
jointLength=accumarray(jointID,weight,[nj 1]);
jointAverage=sparse(jointID,(1:numel(jointID))',weight./jointLength(jointID),nj,numel(jointID));
out.assembly=struct('mesh',M,'parameters',p,'Kelastic',Kelastic, ...
 'thermalLoadMap',thermalMap,'quadrature',quadrature, ...
 'referenceTemperatureGlobal_C',T0(:), ...
 'KringRing',Kfoundation,'KringSoil',coupling.KringSoil,'KsoilSoil',coupling.KsoilSoil, ...
 'Wring',coupling.Wring,'Wsoil',coupling.Wsoil,'Dinterface',coupling.Dinterface, ...
 'prescribedFoundationLoad',foundationLoad, ...
 'prescribedFoundationConstant',foundationConstant, ...
 'prescribedSoilTargetLoad',coupling.soilTargetLoad, ...
 'gapGaussMap',An,'slipGaussMap',At,'gapNodalMap',nodalNormal, ...
 'slipNodalMap',nodalTangent,'gapNodal0_m',repelem(clearance,size(M.jointNodePairs,1)), ...
 'meanGapMap',jointAverage*An,'meanGap0_m',clearance, ...
 'meanSlipMap',jointAverage*At, ...
 'boltNormalMap',Bn,'boltTangentialMap',Bt,'base',base, ...
 'originalGlobalNodeCount',size(M.thermalInterpolation,2), ...
 'mechanicalRingNodeCount',N, ...
 'originalGlobalDisplacementOrder','Interleaved ux(1),uy(1),ux(2),uy(2),...; soil entries only are coupled.', ...
 'soilResidualSign','Internal soil residual Ksr*uring+Kss*usoil; physical force applied to soil is its negative.');
end

function out=evaluateCached(C,state)
% Pure cached algebra for a coupled solver; no mesh/quadrature assembly.
assert(isstruct(C)&&all(isfield(C,{'base','thermalLoadMap','KringSoil','KsoilSoil', ...
 'gapNodalMap','gapNodal0_m'})),'SegmentContactFEM:Cache','Invalid contact FEM cache.');
required={'temperatureGlobal_C','ringDisplacement_m'};
assert(isstruct(state)&&isscalar(state)&&all(isfield(state,required)), ...
 'SegmentContactFEM:State','Cached evaluation needs temperature and ring displacement.');
Ng=C.originalGlobalNodeCount; Nr=C.mechanicalRingNodeCount;
T=state.temperatureGlobal_C(:); Tref=C.referenceTemperatureGlobal_C;
if isfield(state,'referenceTemperatureGlobal_C')
 Tref=state.referenceTemperatureGlobal_C;
 if isscalar(Tref), Tref=repmat(Tref,Ng,1); end
 Tref=Tref(:);
end
assert(numel(T)==Ng&&numel(Tref)==Ng&&all(isfinite(T))&&all(isfinite(Tref)), ...
 'SegmentContactFEM:State','Cached temperature/reference fields have the wrong size.');
ur=displacementVector(state.ringDisplacement_m,Nr); us=zeros(2*Ng,1);
suppliedSoil=isfield(state,'soilDisplacementGlobal_m');
if suppliedSoil
 us=displacementVector(state.soilDisplacementGlobal_m,Ng);
end
extra=zeros(2*Nr,1);
if isfield(state,'additionalRingLoad_N_per_m')
 extra=state.additionalRingLoad_N_per_m(:);
 assert(numel(extra)==2*Nr&&all(isfinite(extra)), ...
  'SegmentContactFEM:State','Additional ring load needs interleaved 2N finite components.');
end
A=C.base;
if suppliedSoil
 foundationLoad=-C.KringSoil*us;
 A.foundationConstant=.5*us'*C.KsoilSoil*us;
 soilTargetLoad=C.KsoilSoil*us;
else
 foundationLoad=C.prescribedFoundationLoad;
 A.foundationConstant=C.prescribedFoundationConstant;
 soilTargetLoad=C.prescribedSoilTargetLoad;
end
A.f=C.thermalLoadMap*(T-Tref)+foundationLoad+extra;
[r,J,energy]=evaluate(ur,A);
soilResidual=C.KringSoil'*ur+soilTargetLoad;
nr=size(C.mesh.jointNodePairs,1); nj=size(C.mesh.jointNodePairs,3);
gap=reshape(C.gapNodal0_m+C.gapNodalMap*ur,nr,nj);
slip=reshape(C.slipNodalMap*ur,nr,nj);
normalStiffness=perJoint(C.parameters.contactNormalStiffness_Pa_per_m,nj);
rotation=zeros(nj,1);
for j=1:nj
 face=C.mesh.jointNodePairs(:,1,j);
 thickness=norm(C.mesh.p(face(end),:)-C.mesh.p(face(1),:));
 rotation(j)=(gap(end,j)-gap(1,j))/thickness;
end
out=struct('ringResidual_N_per_m',r,'ringTangent_N_per_m2',J, ...
 'soilResidual_N_per_m',soilResidual,'ringSoilTangent_N_per_m2',C.KringSoil, ...
 'soilRingTangent_N_per_m2',C.KringSoil','soilSoilTangent_N_per_m2',C.KsoilSoil, ...
 'ringTemperatureTangent_N_per_m_per_C',-C.thermalLoadMap, ...
 'soilReactionForce_global_N_per_m',reshape(-soilResidual,2,[])', ...
 'jointGap_m',gap,'jointSlip_m',slip,'relativeFaceRotation_rad',rotation, ...
 'jointContactPressure_Pa',-min(gap,0).*normalStiffness', ...
 'jointMeanGap_m',C.meanGap0_m+C.meanGapMap*ur, ...
 'jointGapGauss_m',A.g0+A.An*ur, ...
 'jointContactPressureGauss_Pa',-A.contactMaterialK.*min(A.g0+A.An*ur,0), ...
 'boltForce_N_per_m',max(A.preload+A.kb.*(A.Bn*ur),0), ...
 'boltShearForce_N_per_m',A.ks.*(A.Bt*ur), ...
 'potentialEnergy_N',energy,'status','evaluated_cached_contact_and_interface');
if isfield(state,'computeStress')
 assert(isscalar(state.computeStress)&&ismember(state.computeStress,[false true]), ...
  'SegmentContactFEM:State','computeStress must be scalar logical.');
 if state.computeStress
  deltaT=C.mesh.thermalInterpolation*(T-Tref);
  [out.stressGauss_Pa,out.strainGauss]=stressOutput(C.mesh,reshape(ur,2,[])', ...
   deltaT,C.parameters,C.quadrature);
  out.stressComponents={'xx','yy','xy','zz'};
  out.stressCoordinateSystem= ...
   'x,y are original S.pAll coordinates (y is positive-down depth); z is the out-of-plane tunnel axis.';
  out.deltaTemperatureMechanical_C=deltaT;
 end
end
end

function [An,At]=nodalJointOperators(M)
nr=size(M.jointNodePairs,1); nj=size(M.jointNodePairs,3); nd=2*size(M.p,1);
An=sparse(nr*nj,nd); At=An;
for j=1:nj
 normal=[-sin(M.jointAnglesRad(j)),cos(M.jointAnglesRad(j))];
 tangent=[-cos(M.jointAnglesRad(j)),-sin(M.jointAnglesRad(j))];
 for k=1:nr
  dof=udofs(M.jointNodePairs(k,:,j)); row=(j-1)*nr+k;
  An(row,dof)=[-normal,normal]; At(row,dof)=[-tangent,tangent];
 end
end
end

function p=defaults()
p=struct('E_Pa',34.5e9,'poissonRatio',.2,'thermalExpansion_per_C',1e-5, ...
 'referenceTemperatureGlobal_C',20, ...
 'foundationNormalStiffness_Pa_per_m',5e9,'foundationTangentialStiffness_Pa_per_m',1e9, ...
 'foundationDisplacement_m',[],'foundationDisplacementGlobal_m',[], ...
 'contactNormalStiffness_Pa_per_m',1e11,'initialJointGap_m',0, ...
 'boltNormalStiffness_N_per_m2',1e9,'boltPreload_N_per_m',2e4, ...
 'boltShearStiffness_N_per_m2',0,'boltRadialFraction',.75, ...
 'initialDisplacement_m',[],'fixedDofs',[],'fixedValues_m',[], ...
 'solveEquilibrium',true,'relativeTolerance',1e-9, ...
 'absoluteForceTolerance_N_per_m',1e-7,'maximumNewtonIterations',40, ...
 'maximumLineSearchReductions',24,'armijo',1e-4);
end

function validateSettings(p,nd,nj)
assert(isscalar(p.E_Pa)&&isfinite(p.E_Pa)&&p.E_Pa>0&& ...
 isscalar(p.poissonRatio)&&p.poissonRatio>-.99&&p.poissonRatio<.49&& ...
 isscalar(p.thermalExpansion_per_C)&&isfinite(p.thermalExpansion_per_C), ...
 'SegmentContactFEM:Settings','Invalid elastic or thermal-expansion constants.');
nonnegative={'foundationNormalStiffness_Pa_per_m','foundationTangentialStiffness_Pa_per_m'};
for k=1:numel(nonnegative)
 value=p.(nonnegative{k});
 assert(isscalar(value)&&isfinite(value)&&value>=0,'SegmentContactFEM:Settings', ...
  'Foundation stiffness must be a nonnegative scalar.');
end
fields={'contactNormalStiffness_Pa_per_m','initialJointGap_m', ...
 'boltNormalStiffness_N_per_m2','boltPreload_N_per_m','boltShearStiffness_N_per_m2'};
for k=1:numel(fields)
 value=perJoint(p.(fields{k}),nj);
 assert(all(isfinite(value))&&all(value>=0),'SegmentContactFEM:Settings', ...
  'Contact, clearance and bolt parameters must be finite and nonnegative.');
end
position=perJoint(p.boltRadialFraction,nj);
assert(all(isfinite(position))&&all(position>=0&position<=1), ...
 'SegmentContactFEM:Settings','Bolt radial fractions must lie in [0,1].');
assert(isvector(p.fixedDofs)||isempty(p.fixedDofs),'SegmentContactFEM:Settings','Invalid fixed DOFs.');
assert(all(isfinite(p.fixedDofs))&&all(p.fixedDofs==fix(p.fixedDofs))&& ...
 all(p.fixedDofs>=1&p.fixedDofs<=nd)&&numel(unique(p.fixedDofs))==numel(p.fixedDofs), ...
 'SegmentContactFEM:Settings','Fixed DOFs must be unique valid integers.');
assert(isempty(p.fixedValues_m)||(numel(p.fixedValues_m)==numel(p.fixedDofs)&& ...
 all(isfinite(p.fixedValues_m(:)))),'SegmentContactFEM:Settings','Invalid fixed displacement values.');
assert(isempty(p.foundationDisplacement_m)||isempty(p.foundationDisplacementGlobal_m), ...
 'SegmentContactFEM:Settings','Choose mechanical or original global foundation displacements.');
assert(isscalar(p.relativeTolerance)&&p.relativeTolerance>0&& ...
 isscalar(p.absoluteForceTolerance_N_per_m)&&p.absoluteForceTolerance_N_per_m>=0&& ...
 isscalar(p.maximumNewtonIterations)&&p.maximumNewtonIterations>=1&& ...
 p.maximumNewtonIterations==fix(p.maximumNewtonIterations)&& ...
 isscalar(p.maximumLineSearchReductions)&&p.maximumLineSearchReductions>=0&& ...
 p.maximumLineSearchReductions==fix(p.maximumLineSearchReductions)&& ...
 isscalar(p.armijo)&&p.armijo>0&&p.armijo<.5, ...
 'SegmentContactFEM:Settings','Invalid nonlinear-solver settings.');
assert(isscalar(p.solveEquilibrium)&&ismember(p.solveEquilibrium,[false true]), ...
 'SegmentContactFEM:Settings','solveEquilibrium must be scalar logical.');
end

function values=perJoint(values,nj)
assert(isnumeric(values)&&isreal(values)&&(isscalar(values)||numel(values)==nj), ...
 'SegmentContactFEM:Settings','Expected scalar or one parameter per joint.');
values=values(:); if isscalar(values), values=repmat(values,nj,1); end
end

function u=displacementVector(value,N)
if isempty(value), u=zeros(2*N,1); return; end
assert(isnumeric(value)&&isreal(value)&&all(isfinite(value(:)))&& ...
 (isequal(size(value),[N 2])||(isvector(value)&&numel(value)==2*N)), ...
 'SegmentContactFEM:Settings','Displacement must be mechanical N-by-2 or interleaved 2N vector.');
if isequal(size(value),[N 2]), u=reshape(value',[],1); else, u=value(:); end
end

function [K,f,quad,thermalMap]=elasticAssembly(M,dT,p)
nd=2*size(M.p,1); K=sparse(nd,nd); f=zeros(nd,1); ne=size(M.q,1);
mechanicalThermalMap=sparse(nd,size(M.p,1));
D=elasticD(p); eigenstress=p.E_Pa/(1-2*p.poissonRatio)*p.thermalExpansion_per_C*[1;1;0];
quad=struct('B',{cell(ne,4)},'N',{cell(ne,4)});
for e=1:ne
 ids=M.q(e,:); dof=udofs(ids); counter=0;
 for xi=[-1 1]/sqrt(3)
  for eta=[-1 1]/sqrt(3)
   counter=counter+1; [N,d]=qshape(xi,eta); J=d*M.p(ids,:); da=det(J);
   assert(da>0,'SegmentContactFEM:Mesh','Mechanical Q4 Jacobian must be positive.');
   B=bmatrix((J\d)');
   K(dof,dof)=K(dof,dof)+B'*D*B*da;
   f(dof)=f(dof)+B'*eigenstress*(N*dT(ids))*da;
   mechanicalThermalMap(dof,ids)=mechanicalThermalMap(dof,ids)+B'*eigenstress*N*da;
   quad.B{e,counter}=B; quad.N{e,counter}=N;
  end
 end
end
thermalMap=mechanicalThermalMap*M.thermalInterpolation;
end

function [K,f,constant,coupling]=foundationAssembly(M,p)
I=M.interfaceRemapping; nd=2*size(M.p,1); K=sparse(nd,nd); f=zeros(nd,1); constant=0;
nglobal=size(M.thermalInterpolation,2); npoints=2*size(I.mechanicalEdgeNodes,1);
Wring=sparse(2*npoints,nd); Wsoil=sparse(2*npoints,2*nglobal);
Dinterface=sparse(2*npoints,2*npoints); pointIndex=0;
soilTargetLoad=zeros(2*nglobal,1);
assert(isfield(I,'mechanicalEdgeNodes'),'SegmentContactFEM:Mesh','Missing outer mechanical edges.');
mechanical=reshape(displacementVector(p.foundationDisplacement_m,size(M.p,1)),2,[])';
globalField=p.foundationDisplacementGlobal_m;
if ~isempty(globalField)
 assert(I.available&&isequal(size(globalField),[size(M.thermalInterpolation,2) 2])&& ...
  all(isfinite(globalField(:))),'SegmentContactFEM:Settings', ...
  'Original global foundation field needs finite Nglobal-by-2 soil nodal values.');
end
for e=1:size(I.mechanicalEdgeNodes,1)
 ids=I.mechanicalEdgeNodes(e,:); dof=udofs(ids); x=M.p(ids,:);
 tangent=(x(2,:)-x(1,:))'; length=norm(tangent); tangent=tangent/length;
 normal=[tangent(2);-tangent(1)];
 D=p.foundationNormalStiffness_Pa_per_m*(normal*normal')+ ...
  p.foundationTangentialStiffness_Pa_per_m*(tangent*tangent');
 prescribed=mechanical(ids,:);
 if ~isempty(globalField)
  prescribed=squeeze(I.soilEndpointWeights(e,:,:))*globalField(I.soilGlobalNodes(e,:),:);
 end
 for xi=[-1 1]/sqrt(3)
  shape=.5*[1-xi,1+xi]; A=[shape(1)*eye(2),shape(2)*eye(2)];
  target=(shape*prescribed)'; w=length/2;
  K(dof,dof)=K(dof,dof)+A'*D*A*w;
  f(dof)=f(dof)+A'*D*target*w; constant=constant+.5*target'*D*target*w;
  pointIndex=pointIndex+1; rows=2*pointIndex+[-1 0];
  Wring(rows,dof)=A; Dinterface(rows,rows)=w*D;
  if I.available
   soilShape=shape*squeeze(I.soilEndpointWeights(e,:,:));
   As=[soilShape(1)*eye(2),soilShape(2)*eye(2)];
   soilDof=udofs(I.soilGlobalNodes(e,:)); Wsoil(rows,soilDof)=As;
   soilTargetLoad(soilDof)=soilTargetLoad(soilDof)+As'*D*target*w;
  end
 end
end
coupling=struct('KringSoil',-Wring'*Dinterface*Wsoil, ...
 'KsoilSoil',Wsoil'*Dinterface*Wsoil,'Wring',Wring,'Wsoil',Wsoil, ...
 'Dinterface',Dinterface,'soilTargetLoad',soilTargetLoad, ...
 'reactionSign','Physical applied soil reaction is negative soil interface residual.');
end

function [An,At,weight,jointID]=jointGaussOperators(M)
I=M.jointInterface; ng=2*size(I.nodes,1); nd=2*size(M.p,1);
An=sparse(ng,nd); At=An; weight=zeros(ng,1); jointID=weight; row=0;
for e=1:size(I.nodes,1)
 dof=udofs(I.nodes(e,:));
 for xi=[-1 1]/sqrt(3)
  row=row+1; N=.5*[1-xi,1+xi];
  A=[-N(1)*eye(2),-N(2)*eye(2),N(1)*eye(2),N(2)*eye(2)];
  An(row,dof)=I.normalMinusToPlus(e,:)*A;
  At(row,dof)=I.tangentOuterToInner(e,:)*A;
  weight(row)=I.length_m(e)/2; jointID(row)=I.jointID(e);
 end
end
end

function [Bn,Bt]=boltOperators(M,position)
nj=size(M.jointNodePairs,3); nd=2*size(M.p,1); Bn=sparse(nj,nd); Bt=Bn;
for j=1:nj
 minus=M.jointNodePairs(:,1,j); plus=M.jointNodePairs(:,2,j);
 distance=[0;cumsum(hypot(diff(M.p(minus,1)),diff(M.p(minus,2))))];
 target=position(j)*distance(end);
 k=find(distance<=target,1,'last'); k=min(k,numel(distance)-1);
 t=(target-distance(k))/(distance(k+1)-distance(k)); N=[1-t t];
 A=[-N(1)*eye(2),-N(2)*eye(2),N(1)*eye(2),N(2)*eye(2)];
 ids=[minus(k:k+1);plus(k:k+1)]; dof=udofs(ids);
 angle=M.jointAnglesRad(j); normal=[-sin(angle),cos(angle)]; tangent=[-cos(angle),-sin(angle)];
 Bn(j,dof)=normal*A; Bt(j,dof)=tangent*A;
end
end

function [r,J,energy]=evaluate(u,A)
gap=A.g0+A.An*u; overlap=min(gap,0); contact=A.contactK.*overlap;
boltTrial=A.preload+A.kb.*(A.Bn*u); bolt=max(boltTrial,0);
shear=A.ks.*(A.Bt*u);
r=A.K*u-A.f+A.An'*contact+A.Bn'*bolt+A.Bt'*shear;
J=A.K+A.An'*spdiags(A.contactK.*(gap<0),0,numel(gap),numel(gap))*A.An ...
 +A.Bn'*spdiags(A.kb.*(boltTrial>0),0,numel(bolt),numel(bolt))*A.Bn ...
 +A.Bt'*spdiags(A.ks,0,numel(shear),numel(shear))*A.Bt;
active=A.kb>0;
boltEnergy=sum((bolt(active).^2-A.preload(active).^2)./(2*A.kb(active)));
energy=.5*u'*A.K*u-A.f'*u+A.foundationConstant ...
 +.5*sum(A.contactK.*overlap.^2)+boltEnergy+.5*sum(A.ks.*(A.Bt*u).^2);
end

function [gap,slip,rotation]=jointNodalState(M,U,clearance)
nr=size(M.jointNodePairs,1); nj=size(M.jointNodePairs,3); gap=zeros(nr,nj); slip=gap;
rotation=zeros(nj,1);
for j=1:nj
 minus=M.jointNodePairs(:,1,j); plus=M.jointNodePairs(:,2,j); jump=U(plus,:)-U(minus,:);
 angle=M.jointAnglesRad(j); normal=[-sin(angle);cos(angle)]; tangent=[-cos(angle);-sin(angle)];
 opening=jump*normal; gap(:,j)=clearance(j)+opening; slip(:,j)=jump*tangent;
 length=norm(M.p(minus(end),:)-M.p(minus(1),:));
 rotation(j)=(opening(end)-opening(1))/length;
end
end

function [stress,strain]=stressOutput(M,U,dT,p,quad)
ne=size(M.q,1); stress=zeros(ne,4,4); strain=zeros(ne,4,3); D=elasticD(p);
lambda=p.E_Pa*p.poissonRatio/((1+p.poissonRatio)*(1-2*p.poissonRatio));
coefficient=p.E_Pa/(1-2*p.poissonRatio)*p.thermalExpansion_per_C;
for e=1:ne
 ids=M.q(e,:); u=reshape(U(ids,:)',[],1);
 for g=1:4
  epsv=quad.B{e,g}*u; thermal=coefficient*(quad.N{e,g}*dT(ids));
  sigma=D*epsv-thermal*[1;1;0];
  strain(e,g,:)=reshape(epsv,1,1,3);
  stress(e,g,:)=reshape([sigma;lambda*(epsv(1)+epsv(2))-thermal],1,1,4);
 end
end
end

function D=elasticD(p)
E=p.E_Pa; nu=p.poissonRatio;
D=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
end

function dof=udofs(ids)
dof=reshape([2*ids(:)'-1;2*ids(:)'],1,[]);
end

function B=bmatrix(d)
B=zeros(3,2*size(d,1));
for k=1:size(d,1), B(:,2*k-1:2*k)=[d(k,1),0;0,d(k,2);d(k,2),d(k,1)]; end
end

function [N,d]=qshape(xi,eta)
N=.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
d=.25*[-(1-eta),1-eta,1+eta,-(1+eta);-(1-xi),-(1+xi),1+xi,1-xi];
end
