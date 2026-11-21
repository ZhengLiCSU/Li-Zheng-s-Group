function B=build_joint_gap_profiles(C,paths,U)
%BUILD_JOINT_GAP_PROFILES Convex radial contact-gap to thermal-path mapping.
% B=build_joint_gap_profiles(contactCache,thermalPaths) builds B.P mapping
% the contact cache's nr*nj nodal gaps to the EXACT thermal-path sampleXY.
% Rows follow thermal-path order: all nq points of joint 1, then joint 2;
% each path is ordered from outside toward the tunnel. Nodal gap columns
% retain the contact cache's original radial order and independent joint ID.
% Ascending or descending mechanical radii are accepted. Thermal/mechanical
% endpoints and rays must agree; no extrapolation or nearest-point proxy is
% used. B.profileGap0+B.profileGapU*U is a joint-major sample vector.
%
% B=build_joint_gap_profiles(C,paths,U) additionally evaluates nodal/sample
% gap evidence for the supplied ring displacement vector. Without U this
% evidence is explicitly evaluated at U=0, NOT an installed equilibrium.
% Runtime callers can cheaply update it WITHOUT rebuilding this geometry:
% nodal=reshape(B.nodalGap0+B.nodalGapU*U,B.nr,B.nj);
% wholePathOpen=all(nodal>0,1)';
% gapProfile=reshape(B.profileGap0+B.profileGapU*U,B.nq,B.nj)';
%
% A P1 radial gap has a closed continuous cross-section whenever any radial
% node is <=0, including an endpoint missed by midpoint temperature samples.
% Pass wholePathOpen as options.openPathMask to the profile hydraulic kernel
% to block ONLY its new open-slit branch; preserve the base defect and the
% partially open geometric/phase inventories. No positive-gap floor is used.

assert(isstruct(C)&&isscalar(C)&&all(isfield(C,{'mesh','gapNodalMap'})), ...
    'JointGapProfiles:Cache','A contact cache with mesh and gapNodalMap is required.');
M=C.mesh;
assert(isstruct(M)&&isscalar(M)&&all(isfield(M,{'p','center','jointNodePairs'})), ...
    'JointGapProfiles:Mesh','Mechanical mesh requires p, center and jointNodePairs.');
assert(isstruct(paths)&&isscalar(paths)&&all(isfield(paths, ...
    {'nj','nq','sampleXY','radiusOuter_m','radiusInner_m'})), ...
    'JointGapProfiles:Paths','Thermal paths require sampleXY and actual radial limits.');
xy=double(M.p); center=double(M.center(:)'); pairs=double(M.jointNodePairs);
assert(isnumeric(M.p)&&isreal(M.p)&&size(xy,2)==2&&all(isfinite(xy(:))) ...
    &&numel(center)==2&&all(isfinite(center)) ...
    &&size(pairs,2)==2&&all(isfinite(pairs(:)))&&all(pairs(:)==fix(pairs(:))) ...
    &&all(pairs(:)>=1&pairs(:)<=size(xy,1)), ...
    'JointGapProfiles:Mesh','Invalid mechanical coordinates or joint-node inventory.');
nr=size(pairs,1); nj=size(pairs,3); nq=paths.nq; np=nj*nq;
assert(nr>=2&&paths.nj==nj&&isnumeric(nq)&&isscalar(nq)&&isfinite(nq) ...
    &&nq>=1&&nq==fix(nq),'JointGapProfiles:Inventory', ...
    'Thermal and mechanical joint inventories must agree.');
gapMap=sparse(C.gapNodalMap); ng=nr*nj;
assert(isnumeric(C.gapNodalMap)&&isreal(C.gapNodalMap) ...
    &&size(gapMap,1)==ng&&size(gapMap,2)==2*size(xy,1) ...
    &&all(isfinite(nonzeros(gapMap))),'JointGapProfiles:Cache', ...
    'The nodal gap map must use joint-major rows and all ring displacement DOFs.');
if isfield(C,'gapNodal0_m')
    gap0=double(C.gapNodal0_m(:));
elseif isfield(C,'gapNodal0')
    gap0=double(C.gapNodal0(:));
else
    error('JointGapProfiles:Cache','The contact cache requires gapNodal0_m.');
end
assert(isreal(gap0)&&numel(gap0)==ng&&all(isfinite(gap0)), ...
    'JointGapProfiles:Cache','Invalid nodal initial-clearance inventory.');
sample=double(paths.sampleXY); outer=double(paths.radiusOuter_m(:));
inner=double(paths.radiusInner_m(:));
assert(isequal(size(sample),[np 2])&&isreal(sample)&&all(isfinite(sample(:))) ...
    &&numel(outer)==nj&&numel(inner)==nj&&all(isfinite(outer)) ...
    &&all(isfinite(inner))&&all(outer>inner&inner>0), ...
    'JointGapProfiles:Paths','Invalid thermal sample coordinates or limits.');
scale=max(1,max(abs(xy-center),[],'all')); tol=1e-9*scale;
row=zeros(2*np,1); column=row; value=row;
nodalXY=zeros(ng,2); nodalRadius=zeros(nr,nj); sampleRadius=zeros(nq,nj);
for joint=1:nj
    minus=pairs(:,1,joint); plus=pairs(:,2,joint); points=xy(minus,:);
    assert(all(minus~=plus)&&max(abs(points-xy(plus,:)),[],'all')<=tol, ...
        'JointGapProfiles:IndependentFaces', ...
        'Each joint must have coincident faces with independent mechanical nodes.');
    relative=points-center; radius=hypot(relative(:,1),relative(:,2));
    [~,outerNode]=max(radius); direction=relative(outerNode,:)/radius(outerNode);
    if isfield(M,'jointAnglesRad')
        angles=double(M.jointAnglesRad(:));
        assert(numel(angles)==nj&&all(isfinite(angles)), ...
            'JointGapProfiles:Mesh','Invalid joint angle inventory.');
        intended=[cos(angles(joint)),sin(angles(joint))];
        assert(norm(direction-intended,inf)<=tol/max(1,outer(joint)), ...
            'JointGapProfiles:Ray','Mechanical joint angle and node ray disagree.');
    end
    projection=relative*direction';
    assert(all(projection>0)&&max(abs(relative-projection*direction),[],'all')<=tol, ...
        'JointGapProfiles:Ray','Mechanical joint nodes must lie on one outward radial ray.');
    [sortedRadius,order]=sort(projection);
    assert(all(diff(sortedRadius)>tol),'JointGapProfiles:RadialNodes', ...
        'Mechanical radial gap nodes must have distinct radii.');
    assert(abs(sortedRadius(1)-inner(joint))<=tol ...
        &&abs(sortedRadius(end)-outer(joint))<=tol, ...
        'JointGapProfiles:PathBoundary','Thermal and mechanical radial endpoints do not agree.');
    rows=(joint-1)*nq+(1:nq)'; cols=(joint-1)*nr+(1:nr)';
    query=sample(rows,:)-center; rq=query*direction';
    assert(max(abs(query-rq*direction),[],'all')<=tol, ...
        'JointGapProfiles:Ray','Thermal samples do not lie on their mechanical joint ray.');
    assert(all(rq>=sortedRadius(1)-tol&rq<=sortedRadius(end)+tol), ...
        'JointGapProfiles:OutsidePath','Thermal samples are outside the actual mechanical path.');
    assert(all(diff(rq)<0),'JointGapProfiles:SampleOrder', ...
        'Thermal samples must be distinct and ordered outer to inner.');
    nodalXY(cols,:)=points; nodalRadius(:,joint)=projection;
    sampleRadius(:,joint)=rq;
    for k=1:nq
        radiusQuery=rq(k);
        % Values within floating-point tolerance of an endpoint are the
        % exact endpoint; this is not extrapolation of a physical sample.
        radiusQuery=max(sortedRadius(1),min(sortedRadius(end),radiusQuery));
        interval=find(sortedRadius<=radiusQuery,1,'last');
        interval=min(interval,nr-1);
        fraction=(radiusQuery-sortedRadius(interval)) ...
            /(sortedRadius(interval+1)-sortedRadius(interval));
        shape=[1-fraction;fraction]; original=order(interval+(0:1));
        slots=2*(rows(k)-1)+(1:2);
        row(slots)=rows(k); column(slots)=(joint-1)*nr+original;
        value(slots)=shape;
    end
end
P=sparse(row,column,value,np,ng);
assert(all(nonzeros(P)>=0)&&all(abs(full(sum(P,2))-1)<1e-12) ...
    &&max(abs(P*nodalXY-sample),[],'all')<=2*tol, ...
    'JointGapProfiles:Interpolation','Gap interpolation must be convex and reproduce sample coordinates.');
weights=ones(1,nq)/nq;
if isfield(paths,'sampleWeights')&&~isempty(paths.sampleWeights)
    w=paths.sampleWeights;
    assert(isnumeric(w)&&isreal(w)&&isvector(w)&&numel(w)==nq ...
        &&all(isfinite(w(:)))&&all(w(:)>0),'JointGapProfiles:Weights', ...
        'Custom sampleWeights must be positive length fractions.');
    weights=double(w(:)'); weights=weights/sum(weights);
end
if nargin<3 || isempty(U)
    U=zeros(size(gapMap,2),1); evidenceState='U=0 initial-clearance evidence only';
else
    assert(isnumeric(U)&&isreal(U)&&isvector(U)&&numel(U)==size(gapMap,2) ...
        &&all(isfinite(U(:))),'JointGapProfiles:Displacement', ...
        'U must contain every finite ring displacement DOF.');
    U=double(U(:)); evidenceState='supplied ring displacement';
end
nodal=reshape(gap0+gapMap*U,nr,nj); profile0=P*gap0; profileU=P*gapMap;
profile=profile0+profileU*U; gapProfile=reshape(profile,nq,nj)';
wholeOpen=all(nodal>0,1)';
B=struct('P',P,'profileGap0',profile0,'profileGapU',profileU, ...
    'profileGap',profile,'gapProfile_m',gapProfile,'nr',nr,'nj',nj,'nq',nq, ...
    'actualPathLength',outer-inner,'actualPathLength_m',outer-inner, ...
    'sampleWeights',weights,'sampleXY',sample,'nodalXY',nodalXY, ...
    'nodalRadius_m',nodalRadius,'sampleRadius_m',sampleRadius, ...
    'radiusOuter_m',outer,'radiusInner_m',inner, ...
    'nodalGap0',gap0,'nodalGapU',gapMap,'nodalGapEvidence',nodal, ...
    'wholePathOpen',wholeOpen,'evidenceState',evidenceState, ...
    'samplingMissesClosure',all(gapProfile>0,2)&~wholeOpen, ...
    'closureRule','all(reshape(nodalGap0+nodalGapU*U,nr,nj)>0,1)''; block only the new open branch, preserve base defect and inventory', ...
    'sampleOrder','joint-major, outer-to-inner; nodal columns preserve mechanical input order');
B.resistanceIntegration=joint_opening_p1_hydraulics('build',B);
end
