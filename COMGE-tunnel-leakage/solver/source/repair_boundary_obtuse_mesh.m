function report = repair_boundary_obtuse_mesh(sourcePath, outputPath, scope)
%REPAIR_BOUNDARY_OBTUSE_MESH Repair cotangent weights at the lining boundary.
%   REPORT = REPAIR_BOUNDARY_OBTUSE_MESH(SOURCE, OUTPUT) moves only the
%   interior vertices opposite obtuse soil boundary triangles adjacent to
%   the lining. SOURCE is never modified; OUTPUT must not already exist.
%
%   ...(..., 'all') tries the same local projection for every negative
%   boundary edge. It refuses to save if this creates a new negative edge.
%   On the supplied mesh, repairing all 15 edges this way creates 12 new
%   negative interior weights. A connectivity change/remesh is then needed.
%
%   The moved vertices are soil interior points. The lining geometry,
%   soil-lining interface vertices, element connectivity and coupling table
%   remain unchanged. The modified soil triangles change their FEM
%   stiffness and storage geometry, so a fresh sensitivity run is required.

if nargin < 1 || isempty(sourcePath)
    sourcePath = fullfile(fileparts(mfilename('fullpath')), 'example_mesh_verified.mat');
end
if nargin < 2 || isempty(outputPath)
    outputPath = fullfile(fileparts(mfilename('fullpath')), ...
        'example_mesh_repaired_interface.mat');
end
if nargin < 3 || isempty(scope), scope = 'interface'; end
sourcePath = char(sourcePath); outputPath = char(outputPath); scope = char(scope);
assert(any(strcmp(scope, {'interface','all'})), ...
    'Scope must be ''interface'' or ''all''.');
assert(isfile(sourcePath), 'Source mesh does not exist: %s', sourcePath);
sourceFile = java.io.File(sourcePath); outputFile = java.io.File(outputPath);
sourcePath = char(sourceFile.getCanonicalPath());
outputPath = char(outputFile.getCanonicalPath());
assert(~strcmpi(sourcePath, outputPath), 'Output must differ from source mesh.');
assert(~isfile(outputPath), 'Output already exists; no file will be overwritten: %s', outputPath);

S = load(sourcePath);
assert(all(isfield(S, {'pAll','tSoil','soilTunnelCoupleTable'})), ...
    'Mesh lacks pAll, tSoil or soilTunnelCoupleTable.');
xy0 = S.pAll; tri = double(S.tSoil);
assert(size(xy0,2)==2 && size(tri,2)==3 && ...
    all(tri(:)>=1 & tri(:)<=size(xy0,1) & tri(:)==round(tri(:))), ...
    'Invalid soil triangle connectivity.');
interfaceNodes = unique(S.soilTunnelCoupleTable(:,[2 3 6 7]));
interfacePairs = sort(S.soilTunnelCoupleTable(:,[2 3]),2);
[edges, weights0, edgeCount, area20] = edgeWeights(xy0, tri);
assert(all(area20~=0) && all(edgeCount<=2), ...
    'Degenerate or non-manifold input mesh.');
isInterface = ismember(edges, interfacePairs, 'rows');
negative0 = weights0 < -1e-12;
assert(all(edgeCount(negative0)==1), ...
    'Negative interior edge detected; this boundary-only repair is inapplicable.');
if strcmp(scope, 'interface')
    target = find(negative0 & isInterface);
else
    target = find(negative0);
end
assert(~isempty(target), 'No negative boundary edges selected.');

% For a boundary edge AB and its opposite vertex C, cot(angle ACB)>=0
% exactly when C lies outside the circle with diameter AB. Move C only
% along the line from its perpendicular projection onto AB. Keeping AB
% fixed preserves the physical soil-lining and far-field boundaries.
xy = xy0; movedNodes = zeros(numel(target),1); moves = zeros(numel(target),2);
for k = 1:numel(target)
    e = target(k); ends = edges(e,:);
    match = find(sum(ismember(tri, ends),2)==2);
    assert(isscalar(match), 'Selected edge must have exactly one triangle.');
    apex = setdiff(tri(match,:), ends);
    assert(isscalar(apex) && ~ismember(apex,interfaceNodes), ...
        'Opposite vertex is not a soil interior point.');
    assert(~ismember(apex,movedNodes(1:k-1)), ...
        'One vertex opposes multiple selected edges; local repair is ambiguous.');
    A = xy(ends(1),:); B = xy(ends(2),:); C = xy(apex,:);
    edgeVector = B-A; lengthAB = norm(edgeVector); unit = edgeVector/lengthAB;
    along = dot(C-A, unit); foot = A+along*unit;
    normal = C-foot; distance = norm(normal);
    assert(along>0 && along<lengthAB && distance>0, ...
        'Obtuse boundary triangle has unexpected projection geometry.');
    required = sqrt(along*(lengthAB-along));
    margin = 1e-3; % 1 mm margin beyond the right-angle threshold.
    delta = max(0, required+margin-distance)*normal/distance;
    xy(apex,:) = C+delta;
    movedNodes(k) = apex; moves(k,:) = delta;
end

[edges1, weights1, edgeCount1, area21] = edgeWeights(xy, tri);
assert(isequal(edges,edges1) && isequal(edgeCount,edgeCount1), ...
    'Connectivity changed unexpectedly.');
assert(all(area20.*area21>0) && all(abs(area21)>1e-12), ...
    'Repair inverted or collapsed a soil triangle.');
assert(isequal(xy(interfaceNodes,:), xy0(interfaceNodes,:)), ...
    'Interface or lining vertex moved unexpectedly.');
assert(all(weights1(target)>=-1e-12), 'Selected boundary edges remain negative.');
negative1 = weights1 < -1e-12;
newNegative = find(negative1 & ~negative0);
assert(isempty(newNegative), ['Repair created new negative hydraulic edge weights. ' ...
    'No output mesh was written; use a topology-aware remesh.']);
assert(~any(negative1 & isInterface), 'Lining interface has negative edge weights.');

S.pAll = xy;
save(outputPath, '-struct', 'S', '-v7');
report = struct;
report.sourcePath = sourcePath; report.outputPath = outputPath; report.scope = scope;
report.movedNodes = movedNodes; report.displacement_m = moves;
report.maxDisplacement_m = max(vecnorm(moves,2,2));
report.negativeBefore = sum(negative0); report.negativeAfter = sum(negative1);
report.interfaceNegativeBefore = sum(negative0 & isInterface);
report.interfaceNegativeAfter = sum(negative1 & isInterface);
report.newNegativeEdges = numel(newNegative);
report.minimumAreaRatio = min(abs(area21)./abs(area20));
report.remainingNegativeEdges = edges(negative1,:);
end

function [edges, weights, edgeCount, signedArea2] = edgeWeights(xy, tri)
% Match the P1 soil hydraulic stiffness decomposition in TH_2D_segment_v8.
a=xy(tri(:,1),:); b=xy(tri(:,2),:); c=xy(tri(:,3),:);
signedArea2=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2)) ...
    -(c(:,1)-a(:,1)).*(b(:,2)-a(:,2));
assert(all(abs(signedArea2)>1e-12), 'Degenerate soil triangle.');
area=abs(signedArea2)/2;
g1=[b(:,2)-c(:,2),c(:,1)-b(:,1)]./signedArea2;
g2=[c(:,2)-a(:,2),a(:,1)-c(:,1)]./signedArea2;
g3=[a(:,2)-b(:,2),b(:,1)-a(:,1)]./signedArea2;
pairs=[tri(:,[1 2]);tri(:,[2 3]);tri(:,[3 1])];
values=-[area.*sum(g1.*g2,2);area.*sum(g2.*g3,2);area.*sum(g3.*g1,2)];
[edges,~,map]=unique(sort(pairs,2),'rows');
weights=accumarray(map,values,[size(edges,1),1]);
edgeCount=accumarray(map,1,[size(edges,1),1]);
end
