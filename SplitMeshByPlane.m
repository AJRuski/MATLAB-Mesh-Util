function [FacesPositive, VertsPositive, FacesNegative, VertsNegative] = SplitMeshByPlane(Faces, Verts, PlaneNormal, PlanePoint, Options)
%SPLITMESHBYPLANE Cut a watertight mesh with a plane and cap both halves.
%
%   [FacesPositive, VertsPositive, FacesNegative, VertsNegative] = ...
%       SplitMeshByPlane(Faces, Verts, PlaneNormal, PlanePoint)
%   ___ = SplitMeshByPlane(___, Tolerance=1e-9)
%
%   The plane passes through PlanePoint and is perpendicular to
%   PlaneNormal. The "positive" mesh is the part on the side the normal
%   points toward; the "negative" mesh is the part on the other side. Each
%   output mesh is closed off with new triangles lying in the cutting
%   plane, so both are watertight with outward-facing normals.
%
%   How it works:
%     1. Every vertex gets a signed distance to the plane.
%     2. Each edge that crosses the plane gets one new vertex where it
%        crosses. Both triangles that share the edge use the same new
%        vertex, so the pieces stay joined.
%     3. Each triangle that crosses the plane is cut into a smaller
%        triangle on one side and a quadrilateral (two triangles) on the
%        other. The winding of the original triangle is preserved.
%     4. The new vertices form closed loops where the plane meets the
%        surface. The loops are triangulated in 2D with a constrained
%        Delaunay triangulation, which handles non-convex outlines, holes
%        (for example cutting through a tube) and several separate loops.
%        Triangles outside the loops are discarded.
%     5. The cap triangles are added to each piece with the winding that
%        points away from that piece, and unused vertices are dropped.
%
%   The input mesh may be non-convex. A cut through a non-convex mesh can
%   leave a side made of several separate closed pieces; they are returned
%   together in one Faces/Verts pair.
%
%   The mesh must be watertight with outward-facing normals (see
%   IsMeshWatertight and FixFaceNormals); an error is raised otherwise.
%
%   Vertices lying within Tolerance of the plane are treated as being just
%   on the positive side, which avoids zero-length cuts. This can produce
%   very small triangles near such vertices; FixDegenerateTriangles can
%   clean them up if needed.
%
%   Inputs
%       Faces       - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts       - Mx3 matrix of vertex positions.
%       PlaneNormal - 1x3 vector perpendicular to the plane (need not be
%                     unit length).
%       PlanePoint  - 1x3 point on the plane.
%
%   Name-value options
%       Tolerance - Distance from the plane, as a fraction of the mesh
%                   bounding-box diagonal, within which a vertex counts as
%                   on the plane. Default 1e-9.
%
%   Outputs
%       FacesPositive, VertsPositive - Closed mesh on the side PlaneNormal
%                                      points toward.
%       FacesNegative, VertsNegative - Closed mesh on the opposite side.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
    PlaneNormal (1,3) {mustBeNumeric, mustBeReal, mustBeFinite}
    PlanePoint (1,3) {mustBeNumeric, mustBeReal, mustBeFinite}
    Options.Tolerance (1,1) double {mustBePositive, mustBeFinite} = 1e-9
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("SplitMeshByPlane:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

% The plane needs a direction. Normalize it so distances are real lengths.
if norm(PlaneNormal) == 0
    error("SplitMeshByPlane:ZeroNormal", "PlaneNormal must not be the zero vector.");
end
PlaneNormal = PlaneNormal / norm(PlaneNormal);

% The cap winding below assumes outward-facing normals on a closed mesh.
if ~IsMeshWatertight(Faces, Verts)
    error("SplitMeshByPlane:NotWatertight", ...
        "The mesh must be watertight. Use IsMeshWatertight(..., ShowFigure=true) to find the problem.");
end
if ComputeMeshVolume(Faces, Verts) < 0
    error("SplitMeshByPlane:InvertedNormals", ...
        "The mesh normals point inward. Run FixFaceNormals first.");
end

Verts = double(Verts);
VertCount = size(Verts, 1);

%% Step 1: signed distance from each vertex to the plane

Distances = (Verts - PlanePoint) * PlaneNormal';

% Treat vertices that are (almost) on the plane as slightly positive. This
% keeps every distance away from zero, so every crossing edge has a
% well-defined crossing point strictly between its two ends.
BoundingBoxDiagonal = norm(max(Verts, [], 1) - min(Verts, [], 1));
DistanceTolerance = Options.Tolerance * BoundingBoxDiagonal;
Distances(abs(Distances) <= DistanceTolerance) = DistanceTolerance;

IsPositive = Distances > 0;

%% Step 2: one new vertex for every edge that crosses the plane

% Every undirected edge of the mesh, listed once with its smaller vertex
% index first. Sorting the endpoints means both triangles sharing an edge
% find the same entry.
EdgeStart = Faces(:);
EdgeEnd = reshape(Faces(:, [2 3 1]), [], 1);
AllEdges = unique(sort([EdgeStart, EdgeEnd], 2), "rows");

% An edge crosses the plane when its two ends are on different sides.
IsCutEdge = IsPositive(AllEdges(:, 1)) ~= IsPositive(AllEdges(:, 2));
CutEdges = AllEdges(IsCutEdge, :);
CutCount = size(CutEdges, 1);

if CutCount == 0
    error("SplitMeshByPlane:NoIntersection", ...
        "The plane does not cut through the mesh.");
end

% Place each new vertex where its edge crosses the plane. Moving from end
% A toward end B, the distance goes from Da to Db; it is zero after the
% fraction t = Da / (Da - Db) of the way along.
EdgeEndA = CutEdges(:, 1);
EdgeEndB = CutEdges(:, 2);
DistanceA = Distances(EdgeEndA);
DistanceB = Distances(EdgeEndB);
CrossingFraction = DistanceA ./ (DistanceA - DistanceB);
NewVerts = Verts(EdgeEndA, :) + ...
    CrossingFraction .* (Verts(EdgeEndB, :) - Verts(EdgeEndA, :));

% The new vertices are appended after the original ones, so new vertex k
% has index VertCount + k.
AllVerts = [Verts; NewVerts];

%% Step 3: cut the triangles that cross the plane

FacesPositiveClipped = ClipTriangles(Faces, IsPositive, CutEdges, VertCount);
FacesNegativeClipped = ClipTriangles(Faces, ~IsPositive, CutEdges, VertCount);

%% Step 4: find where the plane meets the surface and fill it

% Each triangle that crosses the plane has exactly two cut edges, and the
% new vertices on them are joined by one straight segment of the outline.
PositiveCorners = sum(reshape(IsPositive(Faces), size(Faces)), 2);
CrossesPlane = PositiveCorners == 1 | PositiveCorners == 2;
CrossingFaces = Faces(CrossesPlane, :);

% New-vertex index on each of the three edges (A-B, B-C, C-A) of the
% crossing triangles; 0 where the edge is not cut.
EdgeVertexIds = reshape(LookupCutVertices( ...
    CrossingFaces(:), reshape(CrossingFaces(:, [2 3 1]), [], 1), ...
    CutEdges, VertCount), [], 3);

% Sorting each row in descending order moves the single 0 to the end,
% leaving the two real indices in the first two columns.
EdgeVertexIds = sort(EdgeVertexIds, 2, "descend");
OutlineSegments = unique(sort(EdgeVertexIds(:, 1:2), 2), "rows") - VertCount;

% Describe the cutting plane with a right-handed 2D coordinate system
% (AxisU, AxisV) such that AxisU x AxisV = PlaneNormal. A triangle that
% runs counter-clockwise in (U, V) then has its normal along +PlaneNormal.
[~, WeakestAxis] = min(abs(PlaneNormal));
HelperAxis = zeros(1, 3);
HelperAxis(WeakestAxis) = 1;
AxisU = cross(PlaneNormal, HelperAxis);
AxisU = AxisU / norm(AxisU);
AxisV = cross(PlaneNormal, AxisU);

NewVertsOffset = NewVerts - PlanePoint;
Points2D = [NewVertsOffset * AxisU', NewVertsOffset * AxisV'];

% Two new vertices at exactly the same spot would be merged by the
% triangulation and leave one of them unconnected.
if size(unique(Points2D, "rows"), 1) ~= CutCount
    error("SplitMeshByPlane:CoincidentCutPoints", ...
        "Two cut points coincide. The mesh probably has duplicate vertices or degenerate triangles.");
end

% Triangulate the points with the outline segments forced to be edges, then
% keep only the triangles inside the outline. Holes and separate loops are
% handled by isInterior.
OutlineTriangulation = delaunayTriangulation(Points2D, OutlineSegments);
CapTriangles = OutlineTriangulation.ConnectivityList( ...
    isInterior(OutlineTriangulation), :);

% Make every cap triangle counter-clockwise in (U, V). Twice the signed 2D
% area is the 2D cross product of two edges; negative means clockwise.
CornerA = Points2D(CapTriangles(:, 1), :);
CornerB = Points2D(CapTriangles(:, 2), :);
CornerC = Points2D(CapTriangles(:, 3), :);
TwiceSignedArea = (CornerB(:, 1) - CornerA(:, 1)) .* (CornerC(:, 2) - CornerA(:, 2)) ...
                - (CornerB(:, 2) - CornerA(:, 2)) .* (CornerC(:, 1) - CornerA(:, 1));
IsClockwise = TwiceSignedArea < 0;
CapTriangles(IsClockwise, :) = CapTriangles(IsClockwise, [1 3 2]);

% Convert from numbering within the new vertices to numbering in AllVerts.
CapTriangles = CapTriangles + VertCount;

%% Step 5: close each piece with the cap and tidy up

% The cap normal is +PlaneNormal as built. That is outward for the negative
% piece (which lies below the plane), so it is used as-is there; the
% positive piece lies above the plane and needs the cap reversed.
CapForNegative = CapTriangles;
CapForPositive = CapTriangles(:, [1 3 2]);

[FacesPositive, VertsPositive] = RemoveUnusedVerts( ...
    [FacesPositiveClipped; CapForPositive], AllVerts);
[FacesNegative, VertsNegative] = RemoveUnusedVerts( ...
    [FacesNegativeClipped; CapForNegative], AllVerts);

end


function ClippedFaces = ClipTriangles(Faces, IsKept, CutEdges, VertCount)
%CLIPTRIANGLES Keep the part of every triangle on the kept side of the plane.
%   IsKept is an Mx1 logical saying which original vertices are on the side
%   to keep. Triangles are classified by how many of their corners are kept:
%     3 - kept whole
%     0 - dropped
%     1 - replaced by a small triangle at the kept corner
%     2 - replaced by a quadrilateral (two triangles) around the kept edge
%   The cyclic vertex order of each original triangle is preserved, so the
%   pieces keep the original winding.

% Indexing a column vector with a single-row matrix returns a column, not
% a row, so every lookup below is reshaped back to the shape of the index.
KeptCorners = sum(reshape(IsKept(Faces), size(Faces)), 2);

% Triangles entirely on the kept side are unchanged.
WholeFaces = Faces(KeptCorners == 3, :);

% One kept corner. Rotate each triangle to read [P Q R] with P the kept
% corner. Q and R are then both cut off, and the piece is the small
% triangle [P, crossing on P-Q, crossing on P-R].
OneKept = Faces(KeptCorners == 1, :);
[~, KeptPosition] = max(reshape(IsKept(OneKept), size(OneKept)), [], 2);   % first true in each row
Rotated = RotateRows(OneKept, KeptPosition);
CrossingPQ = LookupCutVertices(Rotated(:, 1), Rotated(:, 2), CutEdges, VertCount);
CrossingPR = LookupCutVertices(Rotated(:, 1), Rotated(:, 3), CutEdges, VertCount);
SmallFaces = [Rotated(:, 1), CrossingPQ, CrossingPR];

% Two kept corners. Rotate each triangle to read [P Q R] with R the cut-off
% corner (so P and Q are kept). The piece is the quadrilateral
% [P, Q, crossing on Q-R, crossing on R-P], split into two triangles.
TwoKept = Faces(KeptCorners == 2, :);
[~, MissingPosition] = min(reshape(IsKept(TwoKept), size(TwoKept)), [], 2);   % first false in each row
Rotated = RotateRows(TwoKept, mod(MissingPosition, 3) + 1);
CrossingQR = LookupCutVertices(Rotated(:, 2), Rotated(:, 3), CutEdges, VertCount);
CrossingRP = LookupCutVertices(Rotated(:, 3), Rotated(:, 1), CutEdges, VertCount);
QuadFirstHalf = [Rotated(:, 1), Rotated(:, 2), CrossingQR];
QuadSecondHalf = [Rotated(:, 1), CrossingQR, CrossingRP];

ClippedFaces = [WholeFaces; SmallFaces; QuadFirstHalf; QuadSecondHalf];

end


function Rotated = RotateRows(Faces, StartPosition)
%ROTATEROWS Cyclically rotate each row so it starts at the given column.
%   A row [A B C] with StartPosition 2 becomes [B C A]. Rotating keeps the
%   triangle's winding.

RowCount = size(Faces, 1);
ColumnIndex = mod(StartPosition - 1 + (0:2), 3) + 1;   % RowCount x 3
RowIndex = repmat((1:RowCount)', 1, 3);
Rotated = Faces(sub2ind(size(Faces), RowIndex, ColumnIndex));

end


function VertexIds = LookupCutVertices(EdgeA, EdgeB, CutEdges, VertCount)
%LOOKUPCUTVERTICES Find the new vertex on each edge A-B that was cut.
%   EdgeA and EdgeB are column vectors of vertex indices. Returns the index
%   (in the extended vertex list) of the crossing vertex on each edge, or 0
%   for edges that are not cut. The edge direction does not matter.

[WasCut, CutRow] = ismember(sort([EdgeA, EdgeB], 2), CutEdges, "rows");
VertexIds = zeros(size(EdgeA));
VertexIds(WasCut) = VertCount + CutRow(WasCut);

end


function [Faces, Verts] = RemoveUnusedVerts(Faces, Verts)
%REMOVEUNUSEDVERTS Drop vertices no face uses and renumber the faces.

UsedVerts = unique(Faces(:));
NewIndex = zeros(size(Verts, 1), 1);
NewIndex(UsedVerts) = 1:numel(UsedVerts);

Verts = Verts(UsedVerts, :);
Faces = reshape(NewIndex(Faces), size(Faces));

end
