function model = build_thermo_interface_ring(S,nq,p,options)
%BUILD_THERMO_INTERFACE_RING Mixed-dimensional verification bridge.
% Uses the existing EIGHT independent Q4 segment meshes, not eight local
% wall springs. Both temperature and displacement live on the split ring.
% Supported, small-strain thermoelastic concrete + bilateral elastic normal
% restraint is an explicit LIMITED verification regime; no contact switch,
% gasket damage, surrounding soil PDE or construction sequence is claimed.
arguments
    S (1,1) struct = make_thm_verification_mesh(16)
    nq (1,1) double {mustBeInteger,mustBePositive} = 4
    p (1,1) struct = joint_thermodynamic_parameters()
    options (1,1) struct = struct()
end
v=struct('referenceTemperature_K',279.15,'normalRestraint_Pam',1e12, ...
    'thermalExpansion_Kinv',1e-5,'heatTransfer_Wm2K',100, ...
    'surfaceHeatTransfer_Wm2K',10,'concreteConductivity_WmK',2, ...
    'concreteCapacity_Jm3K',2400*900,'foundationNormal_Pam',5e9, ...
    'foundationTangential_Pam',1e9,'elasticModulus_Pa',34.5e9,'poissonRatio',0.2);
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(v,names{k}),'ThermoInterface:Option','Unknown option %s.',names{k});
    v.(names{k})=options.(names{k});
end
assert(p.externalTraction_Pa==0,'ThermoInterface:Loading','External traction must be supplied through the host, not the old local spring.');
assert(v.referenceTemperature_K>0 && v.normalRestraint_Pam>=0 && ...
    v.heatTransfer_Wm2K>=0 && v.surfaceHeatTransfer_Wm2K>=0 && ...
    v.concreteConductivity_WmK>0 && v.concreteCapacity_Jm3K>0, ...
    'ThermoInterface:Material','Invalid thermal or restraint coefficients.');
angles=(0:7)'*pi/4;
M=build_segment_mechanical_mesh(S,angles);
C=solve_segment_contact_fem('assemble',M,struct('thermalExpansion_per_C',0, ...
    'boltNormalStiffness_N_per_m2',0,'boltPreload_N_per_m',0, ...
    'foundationNormalStiffness_Pa_per_m',v.foundationNormal_Pam, ...
    'foundationTangentialStiffness_Pa_per_m',v.foundationTangential_Pam, ...
    'E_Pa',v.elasticModulus_Pa,'poissonRatio',v.poissonRatio));
paths=build_joint_thermal_paths(S,angles,nq); gap=build_joint_gap_profiles(C,paths);
n=8*nq; ns=size(M.p,1); weight=repelem(gap.actualPathLength_m/nq,nq);
B=gap.profileGapU;
K=C.Kelastic+C.KringRing+B'*spdiags(weight*v.normalRestraint_Pam,0,n,n)*B;
[~,flag]=chol(K); assert(flag==0,'ThermoInterface:MechanicalDomain','Supported stiffness must be positive definite.');
[Kt,capacity,F]=thermalAssembly(M,v);
strain=sparse(12*size(M.q,1),2*ns);
for e=1:size(M.q,1)
    nodes=M.q(e,:); dofs=reshape([2*nodes-1;2*nodes],1,[]);
    for k=1:4
        rows=(e-1)*12+(k-1)*3+(1:3); strain(rows,dofs)=C.quadrature.B{e,k};
    end
end
minus=reshape(M.jointNodePairs(:,1,:),[],1); plus=reshape(M.jointNodePairs(:,2,:),[],1);
nr=numel(minus); selectMinus=sparse((1:nr)',minus,1,nr,ns);
selectPlus=sparse((1:nr)',plus,1,nr,ns);
Pm=gap.P*selectMinus; Pp=gap.P*selectPlus;
H=spdiags(weight*p.wetAreaFraction*v.heatTransfer_Wm2K,0,n,n)*(Pm+Pp);
[bath,innerBath,outerBath]=bathAssembly(M,v.surfaceHeatTransfer_Wm2K);
faces=zeros(8*(nq+1),2); distance=zeros(size(faces,1),1); bid=zeros(size(distance));
for j=1:8
    cells=(j-1)*nq+(1:nq); rows=(j-1)*(nq+1)+(1:nq+1);
    faces(rows,:)=[[0,cells]',[cells,0]'];
    distance(rows)=gap.actualPathLength_m(j)/nq; distance(rows([1,end]))=distance(rows([1,end]))/2;
    bid(rows(1))=2*j-1; bid(rows(end))=2*j;
end
nf=size(faces,1); D=sparse(n,nf);
for f=1:nf
    if faces(f,1)>0, D(faces(f,1),f)=1; end
    if faces(f,2)>0, D(faces(f,2),f)=-1; end
end
model=struct('weight',weight,'capacity',capacity,'stiffness',K,'thermalStiffness',Kt, ...
    'thermalLoad',F,'gapOperator',B,'displacementPressure',K\(B'*diag(weight*p.wetAreaFraction)), ...
    'displacementTemperature',K\F,'heatExchange',H,'exchangeRowSum',full(sum(H,2)), ...
    'exchangeColumnSum',full(sum(H,1))','bathConductance',bath, ...
    'innerBathConductance',innerBath,'outerBathConductance',outerBath, ...
    'referenceTemperature_K',v.referenceTemperature_K,'conductivityWater',0.6, ...
    'conductivityIce',2.2,'faces',faces,'faceDistance',distance,'boundaryID',bid, ...
    'divergence',D,'mesh',M,'gapProfiles',gap,'options',v,'nq',nq,'strainOperator',strain, ...
    'scope','Eight-block supported thermoelastic ring; bilateral reversible restraint; no soil PDE/contact/damage.');
assert(max(abs(full(sum(Kt,2))))<1e-9 && all(capacity>0),'ThermoInterface:ThermalAssembly', ...
    'Constant temperature must be in the conduction null space; capacity must be positive.');
end

function [K,C,F]=thermalAssembly(M,v)
ns=size(M.p,1); K=sparse(ns,ns); C=zeros(ns,1); F=sparse(2*ns,ns);
thermalStress=v.elasticModulus_Pa/(1-2*v.poissonRatio)*v.thermalExpansion_Kinv*[1;1;0];
for e=1:size(M.q,1)
    ids=M.q(e,:); uids=reshape([2*ids-1;2*ids],1,[]);
    for xi=[-1,1]/sqrt(3)
        for eta=[-1,1]/sqrt(3)
            N=[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)]/4;
            dn=[-(1-eta),1-eta,1+eta,-(1+eta);-(1-xi),-(1+xi),1+xi,1-xi]/4;
            jac=dn*M.p(ids,:); area=det(jac);
            assert(area>0,'ThermoInterface:Mesh','Q4 orientation must be positive.');
            grad=jac\dn; strain=zeros(3,8);
            strain(1,1:2:end)=grad(1,:); strain(2,2:2:end)=grad(2,:);
            strain(3,1:2:end)=grad(2,:); strain(3,2:2:end)=grad(1,:);
            K(ids,ids)=K(ids,ids)+v.concreteConductivity_WmK*(grad'*grad)*area;
            C(ids)=C(ids)+v.concreteCapacity_Jm3K*N'*area;
            F(uids,ids)=F(uids,ids)+strain'*thermalStress*N*area;
        end
    end
end
end

function [bath,inner,outer]=bathAssembly(M,h)
ns=size(M.p,1); inner=zeros(ns,1); outer=inner;
q=M.q; edges=[q(:,[1,2]);q(:,[2,3]);q(:,[3,4]);q(:,[4,1])];
[uniqueEdges,~,group]=unique(sort(edges,2),'rows'); count=accumarray(group,1);
edges=uniqueEdges(count==1,:);
for e=1:size(edges,1)
    ids=edges(e,:); w=norm(M.p(ids(2),:)-M.p(ids(1),:))*h/2;
    if all(M.nodeRingIndex(ids)==1), outer(ids)=outer(ids)+w; end
    if all(M.nodeRingIndex(ids)==max(M.nodeRingIndex)), inner(ids)=inner(ids)+w; end
end
bath=inner+outer;
end
