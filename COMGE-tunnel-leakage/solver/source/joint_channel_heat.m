function [heat,dHeatT,dHeatH,exportCorrection]=joint_channel_heat( ...
    T,soilSelect,paths,flow,baseFlow,headSlope,dFactorT,rhoCw,inletHeatHeadDerivative)
%JOINT_CHANNEL_HEAT Conservative equilibrium water heat exchange in joints.
% flow(j,i) is outward water flow [m2/s] from soil node i into joint j.
% paths.interpolation samples concrete temperature at successive locations
% along each joint, grouped by joint from the inlet toward the tunnel.
% The water equilibrates at each sample. Sensible heat released by cooling
% is a negative contribution to the thermal residual at concrete nodes.
% Joint water/ice storage and latent heat storage are not represented.
%
% baseFlow is flow before multiplying by the temperature factor;
% headSlope is d(flow)/d(h) with that factor held fixed;
% dFactorT is d(factor)/d(T) [nj x N]. The returned derivatives include the
% resulting temperature dependence of water flow and inlet enthalpy.
% The optional ninth input is d(flow*Tsoil)/dh [nj x ns], required for a
% consistent boundary quadrature with non-diagonal head derivatives.
% Without it, the legacy diagonal-node law headSlope*diag(Tsoil) is used.
% rhoCw=rho_water*c_water. exportCorrection changes the existing soil-inlet
% exported water enthalpy to the final joint-sample enthalpy [W/m].
%
% Residual-only call, including the fourth output:
% [heat,~,~,correction]=joint_channel_heat(T,S,paths,flow,[],[],[],rhoCw).
% report=joint_channel_heat('selftest') runs conservative/tangent checks.
if nargin==1&&(ischar(T)||(isstring(T)&&isscalar(T)))&&strcmp(T,'selftest')
    heat=selftest(); dHeatT=[]; dHeatH=[]; exportCorrection=[]; return;
end
assert(nargin==8||nargin==9,'joint_channel_heat requires eight or nine inputs.');
if nargin<9, inletHeatHeadDerivative=[]; end
assert(isnumeric(T)&&isreal(T)&&isvector(T)&&all(isfinite(T(:))), ...
    'T must be a finite real temperature vector.');
T=T(:); N=numel(T); ns=size(soilSelect,1);
assert(isstruct(paths)&&isscalar(paths)&&isfield(paths,'interpolation') ...
    &&isfield(paths,'nj')&&isfield(paths,'nq'),'Invalid joint thermal paths.');
nj=paths.nj; nq=paths.nq; np=nj*nq;
assert(isscalar(nj)&&isfinite(nj)&&nj>=0&&nj==fix(nj) ...
    &&isscalar(nq)&&isfinite(nq)&&nq>=1&&nq==fix(nq), ...
    'Joint/path counts must be nonnegative/positive integers.');
assert(isscalar(rhoCw)&&isreal(rhoCw)&&isfinite(rhoCw)&&rhoCw>0, ...
    'rhoCw must be positive and finite.');
validateMatrix(soilSelect,ns,N,'soilSelect');
validateMatrix(paths.interpolation,np,N,'paths.interpolation');
validateMatrix(flow,nj,ns,'flow');
if ~isempty(inletHeatHeadDerivative)
    validateMatrix(inletHeatHeadDerivative,nj,ns,'inletHeatHeadDerivative');
    inletHeatHeadDerivative=sparse(inletHeatHeadDerivative);
end
soilSelect=sparse(soilSelect); P=sparse(paths.interpolation); flow=sparse(flow);
assert(all(nonzeros(soilSelect)>=0)&&all(abs(full(sum(soilSelect,2))-1)<1e-10), ...
    'soilSelect rows must be convex temperature selectors.');
assert(all(nonzeros(P)>=0)&&all(abs(full(sum(P,2))-1)<1e-10), ...
    'Path interpolation rows must be convex temperature interpolants.');
assert(all(nonzeros(flow)>=0),'Channel flow must be outward/nonnegative.');

firstRows=(1:nq:np)';
R=sparse(firstRows,(1:nj)',ones(nj,1),np,nj);
Repeat=kron(speye(nj),ones(nq,1));
stepDifference=-speye(nq)+sparse((2:nq)',(1:nq-1)',ones(nq-1,1),nq,nq);
D=kron(speye(nj),stepDifference);
tw=soilSelect*T; pathDifference=D*(P*T); q=full(sum(flow,2));
repeatedQ=Repeat*q;
u=R*(flow*tw)+repeatedQ.*pathDifference;
heat=-rhoCw*(P'*u);
exportCorrection=-rhoCw*sum(u);
dHeatT=[]; dHeatH=[];
% Residual callers request output 4 but intentionally omit derivative data.
haveDerivativeInputs=~(isempty(baseFlow)&&isempty(headSlope)&&isempty(dFactorT));
if nargout<2||~haveDerivativeInputs, return; end
validateMatrix(baseFlow,nj,ns,'baseFlow');
validateMatrix(headSlope,nj,ns,'headSlope');
validateMatrix(dFactorT,nj,N,'dFactorT');
baseFlow=sparse(baseFlow); headSlope=sparse(headSlope); dFactorT=sparse(dFactorT);
assert(all(nonzeros(baseFlow)>=0)&&all(nonzeros(headSlope)>=0), ...
    'baseFlow and headSlope must be nonnegative.');
baseInletHeat=baseFlow*tw; baseQ=full(sum(baseFlow,2));
factorAction=R*spdiags(baseInletHeat,0,nj,nj) ...
    +spdiags(pathDifference,0,np,np)*Repeat*spdiags(baseQ,0,nj,nj);
dUT=R*flow*soilSelect+spdiags(repeatedQ,0,np,np)*D*P ...
    +factorAction*dFactorT;
if isempty(inletHeatHeadDerivative)
    inletHeatHeadDerivative=headSlope*spdiags(tw,0,ns,ns);
end
dUH=R*inletHeatHeadDerivative ...
    +spdiags(pathDifference,0,np,np)*Repeat*headSlope;
dHeatT=-rhoCw*(P'*dUT);
dHeatH=-rhoCw*(P'*dUH);
end

function validateMatrix(A,nrow,ncolumn,name)
assert(isnumeric(A)&&isreal(A)&&isequal(size(A),[nrow,ncolumn]) ...
    &&all(isfinite(nonzeros(A))),'Invalid %s dimensions/values.',name);
end

function report=selftest()
% Two joints, four soil temperatures and distinct concrete path nodes.
N=10; ns=4; nj=2; nq=3; rhoCw=1000*4180;
soilSelect=sparse(1:ns,1:ns,1,ns,N);
P=sparse([1 1 2 2 3 4 4 5 5 6], ...
    [5 6 6 7 7 8 9 9 10 10], ...
    [.8 .2 .7 .3 1 .6 .4 .75 .25 1],nj*nq,N);
paths=struct('interpolation',P,'nj',nj,'nq',nq);
T=[8;6;5;9;4;1;-2;7;3;0]; h=[3;2;-1;4];
coefficient=sparse([1 1 2 2],[1 2 3 4],[2e-7 .6e-7 .3e-7 1.4e-7],nj,ns);
entry=[.2 .4 .3 .5;.2 .4 .3 .5];
factorMap=sparse([1 1 1 2 2 2],[5 6 7 8 9 10], ...
    [.10 .05 .02 .03 .07 .08],nj,N);
[flow,baseFlow,headSlope,dFactorT]=testFlows(T,h,coefficient,entry,factorMap);
[heat,JT,JH,correction]=joint_channel_heat(T,soilSelect,paths,flow, ...
    baseFlow,headSlope,dFactorT,rhoCw);
tw=soilSelect*T; profile=reshape(P*T,nq,nj); q=full(sum(flow,2));
expectedCorrection=rhoCw*(sum(q.*profile(end,:)')-sum(flow*tw));
budgetError=abs(sum(heat)-correction);
assert(budgetError<1e-11*max(1,abs(correction)), ...
    'Channel thermal residual does not sum to export correction.');
assert(abs(correction-expectedCorrection)<1e-11*max(1,abs(correction)), ...
    'Channel export correction does not telescope from inlet to outlet.');

uniformT=3.25*ones(N,1);
[uniformHeat,emptyT,emptyH,uniformCorrection]=joint_channel_heat( ...
    uniformT,soilSelect,paths,flow,[],[],[],rhoCw);
assert(norm(uniformHeat,inf)<1e-12&&abs(uniformCorrection)<1e-12, ...
    'Uniform channel/water temperature must produce no sensible heat source.');
assert(isempty(emptyT)&&isempty(emptyH),'Residual-only call built derivatives.');
% Heating one interior concrete node merely moves exchange between points;
% unchanged inlet and final outlet temperatures imply zero total source.
interiorT=uniformT; interiorT(6)=interiorT(6)+4;
[interiorHeat,~,~,interiorCorrection]=joint_channel_heat( ...
    interiorT,soilSelect,paths,flow,[],[],[],rhoCw);
assert(any(interiorHeat>0)&&any(interiorHeat<0) ...
    &&abs(sum(interiorHeat))<1e-11&&abs(interiorCorrection)<1e-11, ...
    'Interior point heat exchange must redistribute heat conservatively.');

directionT=[.2;-.3;.7;.1;-.6;.5;.4;-.2;.8;-.9];
directionH=[.5;-.2;.3;-.7]; epsilon=1e-5;
fp=testFlows(T+epsilon*directionT,h,coefficient,entry,factorMap);
fm=testFlows(T-epsilon*directionT,h,coefficient,entry,factorMap);
hp=joint_channel_heat(T+epsilon*directionT,soilSelect,paths,fp,[],[],[],rhoCw);
hm=joint_channel_heat(T-epsilon*directionT,soilSelect,paths,fm,[],[],[],rhoCw);
fdT=(hp-hm)/(2*epsilon);
errorT=norm(JT*directionT-fdT,inf)/max(1,norm(fdT,inf));
fp=testFlows(T,h+epsilon*directionH,coefficient,entry,factorMap);
fm=testFlows(T,h-epsilon*directionH,coefficient,entry,factorMap);
hp=joint_channel_heat(T,soilSelect,paths,fp,[],[],[],rhoCw);
hm=joint_channel_heat(T,soilSelect,paths,fm,[],[],[],rhoCw);
fdH=(hp-hm)/(2*epsilon);
errorH=norm(JH*directionH-fdH,inf)/max(1,norm(fdH,inf));
assert(errorT<1e-7&&errorH<1e-7, ...
    'Analytic channel heat tangent differs from directional finite differences.');
assert(all(JT(1:ns,:)==0,'all')&&all(JH(1:ns,:)==0,'all'), ...
    'Channel heat must be assigned exclusively to concrete path DOFs.');
report=struct('passed',true,'thermalBudgetError_W_per_m',budgetError, ...
    'temperatureDirectionalRelativeError',errorT, ...
    'headDirectionalRelativeError',errorH, ...
    'uniformHeatNorm_W_per_m',norm(uniformHeat,inf), ...
    'interiorPointTotalHeat_W_per_m',sum(interiorHeat));
disp(report);
end

function [flow,baseFlow,headSlope,dFactorT]=testFlows(T,h,coefficient,entry,A)
% Smooth, nonconstant factor ensures the tangent test exercises TH feedback.
argument=A*T; factor=.6+.25*tanh(argument);
dFactorT=spdiags(.25*(1-tanh(argument).^2),0,numel(factor),numel(factor))*A;
drivingHead=repmat(h(:)',size(coefficient,1),1)-entry;
baseFlow=coefficient.*max(drivingHead,0);
flow=spdiags(factor,0,numel(factor),numel(factor))*baseFlow;
headSlope=spdiags(factor,0,numel(factor),numel(factor)) ...
    *(coefficient.*double(drivingHead>0));
end
