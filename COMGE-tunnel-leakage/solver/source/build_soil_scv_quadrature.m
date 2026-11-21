function Q = build_soil_scv_quadrature(xy,tri,order,subdivisionLevel)
%BUILD_SOIL_SCV_QUADRATURE Positive quadrature of median-dual P1 volumes.
% Each triangular element is divided into six triangles; two belong to
% each vertex's subcontrol volume. P evaluates the ORIGINAL P1 field at
% the quadrature points. Deposit integrates into the actual volume owner,
% rather than replacing a nonlinear phase state by its nodal value.
% order=2: three points per subtriangle; order=5: seven positive points.
if nargin<3, order=5; end
if nargin<4, subdivisionLevel=0; end
validateattributes(subdivisionLevel,{'numeric'},{'scalar','integer','>=',0,'<=',5});
validateattributes(xy,{'double'},{'2d','ncols',2,'finite','real'});
validateattributes(tri,{'double'},{'2d','ncols',3,'integer','positive'});
assert(all(tri(:)<=size(xy,1)),'SCV:connectivity','Invalid triangle index.');
switch order
    case 2
        bary=[2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3];
        weights=ones(3,1)/3;
    case 5
        a=.470142064105115; b=.059715871789770;
        c=.101286507323456; d=.797426985353087;
        bary=[1/3 1/3 1/3;a a b;a b a;b a a;c c d;c d c;d c c];
        weights=[.225;repmat(.132394152788506,3,1);repmat(.125939180544827,3,1)];
        weights=weights/sum(weights);
    otherwise
        error('SCV:order','Supported positive quadrature orders are 2 and 5.');
end
ne=size(tri,1); ns=size(xy,1); ng=numel(weights);
nsub=4^subdivisionLevel; n=6*nsub*ng*ne;
x1=xy(tri(:,1),:); x2=xy(tri(:,2),:); x3=xy(tri(:,3),:);
twiceArea=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) ...
    -(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(twiceArea)/2;
assert(all(area>1e-14),'SCV:degenerate','Degenerate triangle.');
vertex=eye(3); center=ones(1,3)/3;
rows=zeros(n,3); cols=rows; values=rows;
owner=zeros(n,1); element=owner; weight=owner; coordinates=zeros(n,2);
cursor=0;
for i=1:3
    j=mod(i,3)+1; k=mod(i+1,3)+1;
    sub={ [vertex(i,:);(vertex(i,:)+vertex(j,:))/2;center], ...
          [vertex(i,:);center;(vertex(k,:)+vertex(i,:))/2] };
    for level=1:subdivisionLevel
        children=cell(4*numel(sub),1);
        for s=1:numel(sub)
            v=sub{s}; m12=(v(1,:)+v(2,:))/2;
            m23=(v(2,:)+v(3,:))/2; m31=(v(3,:)+v(1,:))/2;
            children(4*s-3:4*s)={ [v(1,:);m12;m31]; ...
                [m12;v(2,:);m23];[m31;m23;v(3,:)];[m12;m23;m31] };
        end
        sub=children;
    end
    for s=1:numel(sub)
        shape=bary*sub{s};
        for g=1:ng
            ids=cursor+(1:ne)'; cursor=cursor+ne;
            rows(ids,:)=repmat(ids,1,3); cols(ids,:)=tri;
            values(ids,:)=repmat(shape(g,:),ne,1);
            owner(ids)=tri(:,i); element(ids)=(1:ne)';
            weight(ids)=area*weights(g)/(6*nsub);
            coordinates(ids,:)=shape(g,1)*x1+shape(g,2)*x2+shape(g,3)*x3;
        end
    end
end
Q.P=sparse(rows(:),cols(:),values(:),n,ns);
Q.deposit=sparse(owner,(1:n)',weight,ns,n);
Q.elementAverage=sparse(element,(1:n)',weight./area(element),ne,n);
Q.coordinates=coordinates; Q.owner=owner; Q.element=element;
Q.weight=weight; Q.area=area; Q.order=order;
Q.subdivisionLevel=subdivisionLevel;
Q.volume=full(Q.deposit*ones(n,1));
Q.note='Median-dual SCV quadrature of continuous P1 fields; positive weights.';
end
