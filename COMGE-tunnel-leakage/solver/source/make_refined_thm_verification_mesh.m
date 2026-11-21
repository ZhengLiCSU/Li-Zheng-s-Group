function S = make_refined_thm_verification_mesh(nAngular,nRadial)
%MAKE_REFINED_THM_VERIFICATION_MESH Controlled annular Q4 refinement family.
% Radius, eight cuts, lining thickness and soil annulus are fixed.
arguments
    nAngular (1,1) double {mustBeInteger,mustBePositive} = 16
    nRadial (1,1) double {mustBeInteger,mustBePositive} = 3
end
assert(nAngular>=16 && mod(nAngular,8)==0,'PhaseSwitch:Mesh','Use at least two elements per segment.');
angles=(0:nAngular-1)'*2*pi/nAngular; center=[0 5];
ring=zeros(nRadial+1,nAngular); xy=zeros(0,2);
for r=1:nRadial+1
    radius=2-(r-1)*(.45/nRadial); ids=size(xy,1)+(1:nAngular); ring(r,:)=ids;
    xy(ids,:)=center+radius*[cos(angles),sin(angles)];
end
soil=zeros(3,nAngular);
for r=1:3
    radius=1+r; ids=size(xy,1)+(1:nAngular); soil(r,:)=ids;
    xy(ids,:)=center+radius*[cos(angles),sin(angles)];
end
q=zeros(nRadial*nAngular,4); cursor=0;
for r=1:nRadial
    for i=1:nAngular
        j=mod(i,nAngular)+1; cursor=cursor+1;
        q(cursor,:)=[ring(r,i),ring(r,j),ring(r+1,j),ring(r+1,i)];
    end
end
tri=zeros(4*nAngular,3); cursor=0;
for r=1:2
    for i=1:nAngular
        j=mod(i,nAngular)+1; cursor=cursor+2;
        tri(cursor-1,:)=[soil(r,i),soil(r+1,i),soil(r+1,j)];
        tri(cursor,:)=[soil(r,i),soil(r+1,j),soil(r,j)];
    end
end
couple=zeros(nAngular,8);
for i=1:nAngular
    j=mod(i,nAngular)+1;
    couple(i,:)=[i,soil(1,i),soil(1,j),i,1,ring(1,i),ring(1,j),i];
end
S=struct('pAll',xy,'tSoil',tri,'qTunnelGlobal',q,'xc',center(1),'yc',center(2), ...
    'soilTunnelCoupleTable',couple,'tunnelInnerRingGlobal',ring(end,:)', ...
    'coordInfo',struct('sectionColumns','X, Z positive down'));
end
