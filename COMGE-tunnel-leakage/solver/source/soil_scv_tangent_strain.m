function M=soil_scv_tangent_strain(v,A,V,kind)
%SOIL_SCV_TANGENT_STRAIN Exact positive-SCV tangent for P1 divergence.
% V maps nodal displacement to constant divergence in each element.
% Group Q.deposit/elementAverage*diag(v)*elementSelect before multiplying
% by V. All physical point derivatives and integration weights are retained.
weighted=reshape(v,A.elements,A.blocks);
if strcmp(kind,'deposit')
    local=(weighted.*A.weight)*A.ownerShape;
    map=sparse(A.ownerRows,A.ownerCols,local(:),A.soilNodes,A.elements);
else
    assert(strcmp(kind,'average'),'SCVTangent:StrainKind', ...
        'Use deposit or average.');
    values=sum(weighted.*A.averageWeight,2);
    map=spdiags(values,0,A.elements,A.elements);
end
M=map*V;
end
