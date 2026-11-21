function [q,ledger,DH,DT]=evaluate_joint_drain(h,T,B,p,waterTableDepth)
%EVALUATE_JOINT_DRAIN Consistent unilateral Robin flux and TH derivatives.
% q is an outward nodal water rate [m2/s] per metre of tunnel axis.
% DH=d(q)/d(h) is consistent, symmetric and positive semidefinite away
% from positive-part switching points. DT includes remote lining thermal
% DOFs through the joint freezing/viscosity factor. At a switching point
% the selected active-set slope is zero.
% ledger.nodeFlow distributes each joint's water to soil nodes; its sum
% over rows equals q. entryHeatHeadDerivative differentiates the integrated
% inlet q*Tsoil term, which cannot use a nodally lumped head derivative.
assert(isnumeric(h)&&isreal(h)&&isvector(h)&&all(isfinite(h(:))), ...
    'evaluate_joint_drain:InvalidHead','h must be a finite real vector.');
assert(isnumeric(T)&&isreal(T)&&isvector(T)&&all(isfinite(T(:))), ...
    'evaluate_joint_drain:InvalidTemperature','T must be a finite real vector.');
h=double(h(:)); T=double(T(:)); ns=numel(h); N=numel(T);
assert(isstruct(B)&&isscalar(B)&&all(isfield(B, ...
    {'jointQuadrature','liningZDrop','soilDepth'})), ...
    'evaluate_joint_drain:InvalidBoundary','Missing quadrature/depth boundary data.');
assert(isstruct(p)&&isscalar(p)&&isfield(p,'leakage'), ...
    'evaluate_joint_drain:InvalidParameters','Missing leakage parameters.');
Q=B.jointQuadrature; phi=Q.phi; nj=Q.nj; nqp=size(phi,1);
assert(size(phi,2)==ns&&all(nonzeros(phi)>=0) ...
    &&all(isfinite(nonzeros(phi)))&&all(abs(full(sum(phi,2))-1)<1e-10), ...
    'evaluate_joint_drain:InvalidQuadrature','Invalid quadrature interpolation.');
assert(numel(Q.jointID)==nqp&&numel(Q.sectorID)==nqp&&numel(Q.weight)==nqp ...
    &&all(Q.jointID>=1&Q.jointID<=nj)&&all(Q.sectorID>=1&Q.sectorID<=Q.nsector) ...
    &&all(isfinite(Q.weight))&&all(Q.weight>0), ...
    'evaluate_joint_drain:InvalidQuadrature','Invalid quadrature inventory/weights.');
jointID=Q.jointID(:); sectorID=Q.sectorID(:); weights=Q.weight(:);
Jmap=sparse((1:nqp)',jointID,1,nqp,nj);
Smap=sparse((1:nqp)',sectorID,1,nqp,Q.nsector);
if isfield(B,'soilSelect')
    soilSelect=sparse(B.soilSelect);
else
    assert(isfield(Q,'soilNodes')&&numel(Q.soilNodes)==ns ...
        &&all(Q.soilNodes>=1&Q.soilNodes<=N), ...
        'evaluate_joint_drain:MissingSoilSelector','Provide B.soilSelect or Q.soilNodes.');
    soilSelect=sparse((1:ns)',Q.soilNodes,1,ns,N);
end
assert(isequal(size(soilSelect),[ns N])&&all(nonzeros(soilSelect)>=0) ...
    &&all(abs(full(sum(soilSelect,2))-1)<1e-10), ...
    'evaluate_joint_drain:InvalidSoilSelector','Invalid soil temperature selector.');
drop=double(B.liningZDrop(:)); depth=double(B.soilDepth(:));
assert(numel(drop)==ns&&numel(depth)==ns&&all(isfinite(drop))&&all(isfinite(depth)), ...
    'evaluate_joint_drain:InvalidDepth','Depth/drop arrays must match soil heads.');
assert(isscalar(waterTableDepth)&&isreal(waterTableDepth)&&isfinite(waterTableDepth), ...
    'evaluate_joint_drain:InvalidWaterTable','Water-table depth must be finite and scalar.');
L=p.leakage;
assert(isfield(L,'interiorHead')&&isscalar(L.interiorHead) ...
    &&isreal(L.interiorHead)&&isfinite(L.interiorHead), ...
    'evaluate_joint_drain:InvalidInteriorHead','Interior pressure head must be scalar.');
cref=jointParameter(L.jointConductance,nj,'jointConductance',0);
entry=jointParameter(L.jointEntryHead,nj,'jointEntryHead',0);
factor=ones(nj,1); dFactorT=sparse(nj,N); freeze=[];
freezeEnabled=isfield(L,'jointFreeze')&&L.jointFreeze.enabled;
if freezeEnabled
    assert(isfield(B,'jointThermalPaths'), ...
        'evaluate_joint_drain:MissingThermalPaths','Enabled freezing requires thermal paths.');
    paths=B.jointThermalPaths; P=paths.interpolation; nq=paths.nq;
    assert(paths.nj==nj&&isequal(size(P),[nj*nq N]), ...
        'evaluate_joint_drain:InvalidThermalPaths','Thermal path dimensions do not match.');
    Tprofile=reshape(P*T,nq,nj)';
    [factor,dFactor,freeze]=joint_freeze_conductance(Tprofile,L.jointFreeze);
    dFactorT=sparse(repelem((1:nj)',nq),(1:nj*nq)', ...
        reshape(dFactor',[],1),nj,nj*nq)*P;
end
ceff=cref.*factor;
driving=h-L.interiorHead+drop;
quadDriving=phi*driving;
headExcess=max(quadDriving-entry(jointID),0);
active=quadDriving>entry(jointID);
baseRate=cref(jointID).*weights.*headExcess;
rate=ceff(jointID).*weights.*headExcess;
slope=ceff(jointID).*weights.*active;
nodeFlow=Jmap'*spdiags(rate,0,nqp,nqp)*phi;
baseFlow=Jmap'*spdiags(baseRate,0,nqp,nqp)*phi;
headSlope=Jmap'*spdiags(slope,0,nqp,nqp)*phi;
q=full(sum(nodeFlow,1))';
DH=phi'*spdiags(slope,0,nqp,nqp)*phi;
DT=baseFlow'*dFactorT;
soilTemperature=soilSelect*T;
entryHeatHeadDerivative=Jmap'*spdiags(slope.*(phi*soilTemperature),0,nqp,nqp)*phi;
reference=depth-waterTableDepth-L.interiorHead+drop;
referenceExcess=max(phi*reference-entry(jointID),0);
ledger=struct('byJoint',full(Jmap'*rate),'bySegment',full(Smap'*rate), ...
    'nodeFlow',nodeFlow,'baseFlow',baseFlow,'headSlope',headSlope, ...
    'dFactorT',dFactorT,'entryHeatHeadDerivative',entryHeatHeadDerivative, ...
    'drivingHead',full(Jmap'*(weights.*quadDriving)), ...
    'referenceHead',full(Jmap'*(weights.*(phi*reference))), ...
    'unrestrictedRate',full(Jmap'*(ceff(jointID).*weights.*referenceExcess)), ...
    'conductance',ceff,'entryHead',entry, ...
    'quadratureHead',quadDriving,'quadratureRate',rate);
if freezeEnabled, ledger.freeze=freeze; end
end

function values=jointParameter(values,nj,name,minimum)
values=double(values(:));
if isempty(values), values=repmat(minimum,nj,1);
elseif isscalar(values), values=repmat(values,nj,1); end
assert(numel(values)==nj&&isreal(values)&&all(isfinite(values))&&all(values>=minimum), ...
    'evaluate_joint_drain:InvalidJointParameter','Invalid %s inventory/value.',name);
end
