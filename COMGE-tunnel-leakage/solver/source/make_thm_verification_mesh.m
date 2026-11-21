function S=make_thm_verification_mesh(n)
%MAKE_THM_VERIFICATION_MESH Small conforming annulus for coupled patch tests.
% Not the research tunnel mesh. The ring is 0.45 m thick, soil radii 2..4 m.
if nargin<1, n=16; end
assert(n>=16&&mod(n,8)==0);
angles=(0:n-1)'*2*pi/n; center=[0 5];
ring=zeros(4,n); xy=zeros(0,2);
for r=1:4
    radius=2-(r-1)*.15; ids=size(xy,1)+(1:n); ring(r,:)=ids;
    xy(ids,:)=center+radius*[cos(angles),sin(angles)];
end
soil=zeros(3,n);
for r=1:3
    radius=1+r; ids=size(xy,1)+(1:n); soil(r,:)=ids;
    xy(ids,:)=center+radius*[cos(angles),sin(angles)];
end
q=zeros(3*n,4); cursor=0;
for r=1:3
    for i=1:n
        j=mod(i,n)+1; cursor=cursor+1;
        q(cursor,:)=[ring(r,i),ring(r,j),ring(r+1,j),ring(r+1,i)];
    end
end
tri=zeros(4*n,3); cursor=0;
for r=1:2
    for i=1:n
        j=mod(i,n)+1; cursor=cursor+2;
        tri(cursor-1,:)=[soil(r,i),soil(r+1,i),soil(r+1,j)];
        tri(cursor,:)=[soil(r,i),soil(r+1,j),soil(r,j)];
    end
end
couple=zeros(n,8);
for i=1:n
    j=mod(i,n)+1;
    couple(i,:)=[i,soil(1,i),soil(1,j),i,1,ring(1,i),ring(1,j),i];
end
S=struct('pAll',xy,'tSoil',tri,'qTunnelGlobal',q,'xc',center(1),'yc',center(2), ...
    'soilTunnelCoupleTable',couple,'tunnelInnerRingGlobal',ring(end,:)', ...
    'coordInfo',struct('sectionColumns','X, Z positive down'));
end
