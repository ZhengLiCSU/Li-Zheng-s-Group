function A=build_scv_tangent_assembly(Q,tri)
%BUILD_SCV_TANGENT_ASSEMBLY Fixed layouts for exact positive-SCV tangents.
% Q stores quadrature block by block, each block containing all elements.
% Grouping the same weighted products into element-local matrices avoids
% repeatedly constructing diagonal sparse matrices. No quadrature point,
% shape function, physical derivative or integration weight is removed.
ne=size(tri,1); ns=size(Q.P,2); nqp=size(Q.P,1);
assert(mod(nqp,ne)==0,'SCVTangent:Layout','Expected equal quadrature blocks per element.');
blocks=nqp/ne; first=(1:ne:nqp)';
shape=full(Q.P(first,tri(1,:)));
assert(all(abs(sum(shape,2)-1)<1e-13),'SCVTangent:Shape','Invalid P1 block layout.');
owners=Q.owner(first); ownerLocal=zeros(blocks,1);
for i=1:3
    ownerLocal(owners==tri(1,i))=i;
end
assert(all(ownerLocal>0),'SCVTangent:Owner','Unknown volume owner.');
assert(isequal(reshape(Q.owner,ne,blocks),tri(:,ownerLocal)), ...
    'SCVTangent:Layout','Quadrature owners are not grouped by element.');
liftShape=zeros(blocks,9); rows=zeros(ne,9); cols=rows;
ownerShape=zeros(blocks,3);
for i=1:3
    ownerShape(ownerLocal==i,i)=1;
    for j=1:3
        k=(i-1)*3+j;
        liftShape(ownerLocal==i,k)=shape(ownerLocal==i,j);
        rows(:,k)=tri(:,i); cols(:,k)=tri(:,j);
    end
end
A=struct('elements',ne,'soilNodes',ns,'blocks',blocks, ...
    'weight',reshape(Q.weight,ne,blocks), ...
    'averageWeight',reshape(Q.weight./Q.area(Q.element),ne,blocks), ...
    'liftShape',liftShape,'averageShape',shape,'ownerShape',ownerShape, ...
    'liftRows',rows(:),'liftCols',cols(:), ...
    'averageRows',repmat((1:ne)',3,1),'averageCols',tri(:), ...
    'ownerRows',tri(:),'ownerCols',repmat((1:ne)',3,1));
end
