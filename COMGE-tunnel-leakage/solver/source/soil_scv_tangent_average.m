function M=soil_scv_tangent_average(v,A)
%SOIL_SCV_TANGENT_AVERAGE Same Q.elementAverage*diag(v)*Q.P.
local=(reshape(v,A.elements,A.blocks).*A.averageWeight)*A.averageShape;
M=sparse(A.averageRows,A.averageCols,local(:),A.elements,A.soilNodes);
end
