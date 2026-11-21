function Q=build_joint_drain_quadrature(S,soilNodes,leakage,segmentBoundaryAnglesRad)
%BUILD_JOINT_DRAIN_QUADRATURE Consistent finite-band drainage on P1 edges.
% Each band carries one integrated joint conductance per metre of tunnel
% axis. Two-point Gauss integration preserves the consistent Robin matrix,
% including off-diagonal head derivatives, instead of nodal mass lumping.
% The band width locates the soil drainage footprint; it is not an opening.
% Cuts at band and segment angles use exact intersections with polygonal
% soil-interface edges. phi interpolates soil fields at the Gauss points.
assert(isstruct(S)&&isscalar(S)&&all(isfield(S, ...
    {'pAll','xc','yc','soilTunnelCoupleTable'})), ...
    'build_joint_drain_quadrature:InvalidMesh','Missing mesh interface data.');
xy=double(S.pAll); center=double([S.xc S.yc]);
assert(size(xy,2)==2&&all(isfinite(xy(:)))&&numel(center)==2 ...
    &&all(isfinite(center)), ...
    'build_joint_drain_quadrature:InvalidCoordinates','Invalid mesh coordinates.');
soilNodes=double(soilNodes(:)); ns=numel(soilNodes);
assert(all(isfinite(soilNodes))&&all(soilNodes==fix(soilNodes)) ...
    &&all(soilNodes>=1&soilNodes<=size(xy,1)) ...
    &&numel(unique(soilNodes))==ns, ...
    'build_joint_drain_quadrature:InvalidSoilNodes','Invalid soil node inventory.');
assert(isstruct(leakage)&&isscalar(leakage)&& ...
    all(isfield(leakage,{'jointAnglesRad','jointWidth_m'})), ...
    'build_joint_drain_quadrature:InvalidLeakage','Missing joint angle/width data.');
angles=double(leakage.jointAnglesRad(:)); nj=numel(angles);
assert(all(isfinite(angles))&&isreal(angles), ...
    'build_joint_drain_quadrature:InvalidAngles','Joint angles must be finite and real.');
widths=double(leakage.jointWidth_m(:));
if isscalar(widths), widths=repmat(widths,nj,1); end
assert(numel(widths)==nj&&isreal(widths)&&all(isfinite(widths))&&all(widths>0), ...
    'build_joint_drain_quadrature:InvalidWidth','Supply one positive width per joint.');
bounds=double(segmentBoundaryAnglesRad(:)); nsector=numel(bounds);
assert(nsector>=1&&all(isfinite(bounds))&&isreal(bounds), ...
    'build_joint_drain_quadrature:InvalidSectors','Segment boundaries must be finite angles.');
origin=bounds(1); relativeBounds=mod(bounds-origin,2*pi);
assert(relativeBounds(1)<1e-12&&all(diff(relativeBounds)>1e-10) ...
    &&relativeBounds(end)<2*pi-1e-10, ...
    'build_joint_drain_quadrature:SectorOrder', ...
    'Segment boundaries must follow their circular order, without duplicates.');
table=S.soilTunnelCoupleTable;
assert(size(table,2)>=3, ...
    'build_joint_drain_quadrature:InvalidInterface','Interface table lacks edge columns.');
edges=double(table(:,2:3));
assert(all(isfinite(edges(:)))&&all(edges(:)==fix(edges(:))) ...
    &&all(edges(:)>=1&edges(:)<=size(xy,1)), ...
    'build_joint_drain_quadrature:InvalidInterface','Invalid interface node IDs.');
map=zeros(size(xy,1),1); map(soilNodes)=1:ns;
assert(all(map(edges(:))>0), ...
    'build_joint_drain_quadrature:NonSoilInterface','Interface must use soil-side nodes.');
assert(size(unique(sort(edges,2),'rows'),1)==size(edges,1), ...
    'build_joint_drain_quadrature:RepeatedInterface','Repeated interface edges.');
if nj==0
    Q=struct('phi',sparse(0,ns),'jointID',zeros(0,1), ...
        'sectorID',zeros(0,1),'weight',zeros(0,1), ...
        'supportLength',zeros(0,1),'sampleXY',zeros(0,2), ...
        'soilNodes',soilNodes,'nj',0,'nsector',nsector);
    return;
end
assert(~isempty(edges), ...
    'build_joint_drain_quadrature:EmptyInterface','Joint drainage needs an interface.');
boundaryNodes=unique(edges(:));
radius=mean(hypot(xy(boundaryNodes,1)-center(1),xy(boundaryNodes,2)-center(2)));
assert(isfinite(radius)&&radius>0&&all(widths/(2*radius)<pi/2), ...
    'build_joint_drain_quadrature:ExcessiveWidth','Each band must span less than half a ring.');
for first=1:nj
    for second=first+1:nj
        separation=abs(mod(angles(first)-angles(second)+pi,2*pi)-pi);
        assert(radius*separation>=(widths(first)+widths(second))/2-1e-10, ...
            'build_joint_drain_quadrature:OverlappingBands', ...
            'Overlapping bands require a single combined conductance.');
    end
end
row=zeros(0,1); column=row; values=row;
jointID=row; sectorID=row; rawWeight=row; sampleXY=zeros(0,2);
gauss=[.5-.5/sqrt(3);.5+.5/sqrt(3)];
for joint=1:nj
    halfAngle=widths(joint)/(2*radius);
    for edge=1:size(edges,1)
        nodes=edges(edge,:); a=xy(nodes(1),:)-center; b=xy(nodes(2),:)-center;
        theta=[atan2(a(2),a(1)),atan2(b(2),b(1))];
        theta(2)=theta(1)+mod(theta(2)-theta(1)+pi,2*pi)-pi;
        assert(abs(diff(theta))>1e-12&&abs(diff(theta))<pi-1e-10, ...
            'build_joint_drain_quadrature:DegenerateAngularEdge', ...
            'Interface edges must have a nondegenerate angular span below pi.');
        near=angles(joint)+2*pi*round((mean(theta)-angles(joint))/(2*pi));
        low=max(min(theta),near-halfAngle); high=min(max(theta),near+halfAngle);
        if high<=low+1e-14, continue; end
        copies=reshape(mod(bounds,2*pi)+2*pi*(-2:2),[],1);
        crossings=copies(copies>low+1e-12&copies<high-1e-12);
        cuts=sort([low;crossings;high]);
        v=b-a; edgeLength=norm(v);
        for piece=1:numel(cuts)-1
            u1=rayFraction(a,v,cuts(piece)); u2=rayFraction(a,v,cuts(piece+1));
            uLow=min(u1,u2); uHigh=max(u1,u2);
            fractions=uLow+(uHigh-uLow)*gauss;
            positions=xy(nodes(1),:)+fractions*v;
            mid=mean(positions,1)-center;
            angle=mod(atan2(mid(2),mid(1))-origin,2*pi);
            sector=find(relativeBounds<=angle+1e-12,1,'last');
            if isempty(sector), sector=nsector; end
            indices=numel(jointID)+(1:2)';
            jointID(indices,1)=joint; sectorID(indices,1)=sector; %#ok<AGROW>
            rawWeight(indices,1)=edgeLength*(uHigh-uLow)/2; %#ok<AGROW>
            sampleXY(indices,:)=positions; %#ok<AGROW>
            row=[row;indices;indices]; %#ok<AGROW>
            column=[column;repmat(map(nodes(1)),2,1);repmat(map(nodes(2)),2,1)]; %#ok<AGROW>
            values=[values;1-fractions;fractions]; %#ok<AGROW>
        end
    end
end
supportLength=accumarray(jointID,rawWeight,[nj 1]);
assert(all(supportLength>0), ...
    'build_joint_drain_quadrature:UnsupportedJoint','A joint band misses the interface.');
weight=rawWeight./supportLength(jointID);
phi=sparse(row,column,values,numel(jointID),ns);
assert(all(abs(full(sum(phi,2))-1)<1e-12)&&all(nonzeros(phi)>=0) ...
    &&all(abs(accumarray(jointID,weight,[nj 1])-1)<1e-12), ...
    'build_joint_drain_quadrature:InvalidWeights','Quadrature normalization failed.');
Q=struct('phi',phi,'jointID',jointID,'sectorID',sectorID, ...
    'weight',weight,'supportLength',supportLength,'sampleXY',sampleXY, ...
    'soilNodes',soilNodes,'nj',nj,'nsector',nsector);
end

function fraction=rayFraction(a,v,theta)
direction=[cos(theta),sin(theta)];
denominator=direction(1)*v(2)-direction(2)*v(1);
assert(abs(denominator)>1e-14*max(1,norm(v)), ...
    'build_joint_drain_quadrature:RadialEdge','Band cut is parallel to an interface edge.');
fraction=(a(1)*direction(2)-a(2)*direction(1))/denominator;
assert(fraction>=-1e-9&&fraction<=1+1e-9, ...
    'build_joint_drain_quadrature:InvalidIntersection','Angular cut lies outside its edge.');
fraction=max(0,min(1,fraction));
end
