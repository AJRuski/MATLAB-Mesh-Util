function [Faces, Verts] = JoinMeshes(FacesList, VertsList, Options)
%JOINMESHES Join overlapping meshes (open or closed) into one outer skin.
%
%   [Faces, Verts] = JoinMeshes({FacesA, FacesB, ...}, {VertsA, VertsB, ...})
%   [Faces, Verts] = JoinMeshes(___, Tolerance=1e-9)
%
%   The meshes are cut along the curves where they cross and welded where
%   they share vertices. The result is the outer skin of the joined
%   surfaces: every part of a surface that borders the outside is kept, and
%   every part hidden inside (for example the end of a wing that pokes into
%   a fuselage) is removed. If the parts together enclose a volume, the
%   result is a single watertight mesh with outward-facing normals.
%
%   The input meshes do not have to be watertight. A fuselage with open
%   ends and separate end caps, or a wing made of open surfaces, can be
%   joined as long as the pieces close up where they meet. The winding
%   (normal direction) of the inputs does not matter either; the output
%   normals are rebuilt to point outward.
%
%   How it works:
%     1. Triangles from different meshes whose bounding boxes overlap are
%        tested in pairs. Each pair that crosses meets along one segment.
%     2. Segment end points and vertices within Tolerance of each other
%        are merged, so triangles that touch the same point share one
%        vertex. This also welds meshes that simply touch along an edge.
%     3. Every triangle that has segments in it is retriangulated in 2D
%        with a constrained Delaunay triangulation that forces the
%        segments to become edges, keeping the triangle's winding.
%     4. The cut surfaces divide space into regions. Around every edge the
%        triangles meeting there are sorted by angle, and the two sides of
%        neighboring triangles that face the same wedge are linked. Linked
%        sides form one region.
%     5. The ray from far outside toward the mesh first hits the outside
%        of the surface. The region on that side of the first triangle hit
%        is the outside region. Triangles with a side on the outside
%        region are kept, flipped if needed so that side is the front.
%
%   Limitations:
%     - Faces of different meshes that lie in the same plane and overlap
%       (for example two boxes with a flush side) are not resolved unless
%       their triangles are identical. Moving one mesh a tiny amount avoids
%       the problem.
%     - Open edges must meet another surface exactly. A mesh whose open edge
%       ends partway along a triangle edge of another mesh (a T-junction)
%       is not stitched, and the result will have a hole there.
%     - A group of triangles that is not connected to the rest, such as a
%       whole mesh floating inside another, is kept as its own piece and a
%       warning is issued.
%     - A surface that is exposed on both sides (a loose sheet) is kept as
%       a sheet, so the result will not be watertight.
%     - Triangles of one mesh are not tested against each other, so a
%       self-intersecting mesh stays self-intersecting.
%     - Every triangle is compared with every other, so meshes with
%       hundreds of thousands of triangles will be slow.
%
%   Inputs
%       FacesList - Cell array of Nx3 triangle index matrices, one per mesh.
%       VertsList - Cell array of Mx3 position matrices, one per mesh, in
%                   the same order as FacesList.
%
%   Name-value options
%       Tolerance - Distance, as a fraction of the combined bounding-box
%                   diagonal, within which two points count as the same
%                   point or a point counts as lying on a plane or an edge.
%                   Default 1e-9.
%
%   Outputs
%       Faces - Kx3 triangle vertex indices of the joined mesh.
%       Verts - Px3 vertex positions of the joined mesh.

arguments
    FacesList (1,:) cell {mustBeNonempty}
    VertsList (1,:) cell {mustBeNonempty}
    Options.Tolerance (1,1) double {mustBePositive, mustBeFinite} = 1e-9
end

if numel(FacesList) ~= numel(VertsList)
    error("JoinMeshes:CountMismatch", ...
        "FacesList has %d meshes but VertsList has %d.", numel(FacesList), numel(VertsList));
end

%% Step 0: put all the meshes into one list

MeshCount = numel(FacesList);
FaceBlocks = cell(MeshCount, 1);
VertBlocks = cell(MeshCount, 1);
MeshOfFaceBlocks = cell(MeshCount, 1);
VertOffset = 0;

for MeshIndex = 1:MeshCount
    MeshFaces = FacesList{MeshIndex};
    MeshVerts = VertsList{MeshIndex};
    validateattributes(MeshFaces, {'numeric'}, {'2d', 'ncols', 3, 'nonempty', 'real', 'integer', 'positive'}, ...
        "JoinMeshes", sprintf("FacesList{%d}", MeshIndex));
    validateattributes(MeshVerts, {'numeric'}, {'2d', 'ncols', 3, 'nonempty', 'real', 'finite'}, ...
        "JoinMeshes", sprintf("VertsList{%d}", MeshIndex));

    % Every index must point at an existing vertex.
    if any(MeshFaces(:) > size(MeshVerts, 1))
        error("JoinMeshes:InvalidFaceIndex", ...
            "FacesList{%d} references a vertex index larger than the %d vertices provided.", ...
            MeshIndex, size(MeshVerts, 1));
    end

    % Vertex numbers of later meshes are shifted past the earlier ones.
    FaceBlocks{MeshIndex} = double(MeshFaces) + VertOffset;
    VertBlocks{MeshIndex} = double(MeshVerts);
    MeshOfFaceBlocks{MeshIndex} = repmat(MeshIndex, size(MeshFaces, 1), 1);
    VertOffset = VertOffset + size(MeshVerts, 1);
end

AllFaces = vertcat(FaceBlocks{:});
AllVerts = vertcat(VertBlocks{:});
MeshOfFace = vertcat(MeshOfFaceBlocks{:});
VertTotal = size(AllVerts, 1);

% Work with absolute distances from here on.
BoundingBoxDiagonal = norm(max(AllVerts, [], 1) - min(AllVerts, [], 1));
Tolerance = Options.Tolerance * BoundingBoxDiagonal;

%% Step 1: find where the surfaces of different meshes cross

% Each row of the pair list is two crossing triangles (numbers into
% AllFaces) with the end points of the segment where they meet.
[PairFaceI, PairFaceJ, SegmentStart, SegmentEnd] = ...
    FindIntersectionSegments(AllFaces, AllVerts, MeshOfFace, Tolerance);
SegmentCount = numel(PairFaceI);

%% Step 2: merge the segment end points with each other and with the vertices

% One shared list holds every vertex, then every segment end point. Points
% within Tolerance of each other collapse to one entry, and MergeMap says
% which merged vertex each original entry became. The segment end points
% are what make neighboring triangles, and the meshes, agree on a vertex.
AllPoints = [AllVerts; SegmentStart; SegmentEnd];
[MergedVerts, ~, MergeMap] = uniquetol(AllPoints, Tolerance, ...
    "ByRows", true, "DataScale", 1);
MergeMap = MergeMap(:);

MergedFaces = reshape(MergeMap(AllFaces), size(AllFaces));
SegmentStartId = MergeMap(VertTotal + (1:SegmentCount)');
SegmentEndId = MergeMap(VertTotal + SegmentCount + (1:SegmentCount)');
Segments = [SegmentStartId, SegmentEndId];

%% Step 3: cut the triangles along the segments

% Each segment lies in both of the triangles that produced it.
CutFaces = SplitIntersectedTriangles(MergedFaces, MergedVerts, ...
    [PairFaceI; PairFaceJ], [Segments; Segments], Tolerance);

% Tidy up: merging nearby vertices can collapse tiny triangles, and the
% region search below needs every triangle to have a real direction.
CutFaces = RemoveRepeatedIndexFaces(CutFaces);
CutFaces = RemoveZeroAreaFaces(CutFaces, MergedVerts, Tolerance ^ 2);
CutFaces = RemoveDuplicateFaces(CutFaces);

if isempty(CutFaces)
    error("JoinMeshes:NoFaces", "No triangles are left after cutting the meshes.");
end

%% Step 4: group the triangle sides into regions of space

% SideRegion(2*k-1) is the region in front of triangle k (where its normal
% points) and SideRegion(2*k) the region behind it. FaceComponent says which
% connected group of triangles each triangle belongs to.
[SideRegion, FaceComponent] = FindSideRegions(CutFaces, MergedVerts);

%% Step 5: keep the triangles that border the outside

KeepFace = false(size(CutFaces, 1), 1);
FlipFace = false(size(CutFaces, 1), 1);
ComponentCount = max(FaceComponent);

for Component = 1:ComponentCount
    ComponentFaces = find(FaceComponent == Component);

    % The region on the outside of this group of triangles.
    ExteriorSide = FindExteriorSide(CutFaces(ComponentFaces, :), MergedVerts, ComponentFaces);
    ExteriorRegion = SideRegion(ExteriorSide);

    % A triangle with the outside region on its front is kept as is; with
    % it on the back it is kept and flipped so the front faces outward.
    FrontIsOutside = SideRegion(2 * ComponentFaces - 1) == ExteriorRegion;
    BackIsOutside = SideRegion(2 * ComponentFaces) == ExteriorRegion;
    KeepFace(ComponentFaces) = FrontIsOutside | BackIsOutside;
    FlipFace(ComponentFaces) = BackIsOutside & ~FrontIsOutside;
end

Skin = CutFaces(KeepFace, :);
IsFlipped = FlipFace(KeepFace);
Skin(IsFlipped, :) = Skin(IsFlipped, [1 3 2]);

[Faces, Verts] = RemoveUnusedVerts(Skin, MergedVerts);

if ComponentCount > 1
    warning("JoinMeshes:SeparatePieces", ...
        "The meshes form %d separate pieces. A piece lying completely inside another is not removed.", ...
        ComponentCount);
end
if ~IsMeshWatertight(Faces, Verts)
    warning("JoinMeshes:ResultNotWatertight", ...
        "The joined mesh is not watertight. Inspect it with IsMeshWatertight(..., ShowFigure=true).");
end

end


function [PairFaceI, PairFaceJ, SegmentStart, SegmentEnd] = FindIntersectionSegments(Faces, Verts, MeshOfFace, Tolerance)
%FINDINTERSECTIONSEGMENTS Find every pair of crossing triangles from different meshes.
%   Returns Kx1 triangle numbers (PairFaceI < PairFaceJ) and Kx3 start and
%   end points of the segment shared by each pair.

[BoxMin, BoxMax] = ComputeFaceBoundingBoxes(Faces, Verts);
FaceCount = size(Faces, 1);

% One entry per triangle, each holding rows of
% [FaceI FaceJ StartXYZ EndXYZ] for the pairs found with a later triangle.
FoundPerFace = repmat({zeros(0, 8)}, FaceCount, 1);

for FaceI = 1:FaceCount - 1

    % Each pair is tested once, so only later triangles are compared.
    LaterFaces = (FaceI + 1:FaceCount)';

    % Cheap rejection first: only triangles of another mesh whose bounding
    % boxes overlap (with a little slack) can possibly cross.
    BoxesOverlap = all(BoxMin(LaterFaces, :) <= BoxMax(FaceI, :) + Tolerance, 2) & ...
                   all(BoxMax(LaterFaces, :) >= BoxMin(FaceI, :) - Tolerance, 2);
    IsOtherMesh = MeshOfFace(LaterFaces) ~= MeshOfFace(FaceI);
    Candidates = LaterFaces(BoxesOverlap & IsOtherMesh);

    TriangleI = Verts(Faces(FaceI, :), :);
    Found = zeros(0, 8);
    for FaceJ = Candidates'
        TriangleJ = Verts(Faces(FaceJ, :), :);
        [DoesCross, Start, Finish] = IntersectTriangles(TriangleI, TriangleJ, Tolerance);
        if DoesCross
            Found(end + 1, :) = [FaceI, FaceJ, Start, Finish]; %#ok<AGROW>
        end
    end
    FoundPerFace{FaceI} = Found;
end

Found = vertcat(FoundPerFace{:});
PairFaceI = Found(:, 1);
PairFaceJ = Found(:, 2);
SegmentStart = Found(:, 3:5);
SegmentEnd = Found(:, 6:8);

end


function [BoxMin, BoxMax] = ComputeFaceBoundingBoxes(Faces, Verts)
%COMPUTEFACEBOUNDINGBOXES Per-triangle minimum and maximum corner (Nx3 each).

CornerA = Verts(Faces(:, 1), :);
CornerB = Verts(Faces(:, 2), :);
CornerC = Verts(Faces(:, 3), :);
BoxMin = min(min(CornerA, CornerB), CornerC);
BoxMax = max(max(CornerA, CornerB), CornerC);

end


function [DoesCross, Start, Finish] = IntersectTriangles(TriangleA, TriangleB, Tolerance)
%INTERSECTTRIANGLES Segment where two triangles cross.
%   TriangleA and TriangleB are 3x3 (one corner per row). DoesCross is false
%   when the triangles miss, only touch at a point, or are coplanar.
%
%   Each triangle meets the other's plane along a short stretch of the line
%   where the two planes cross. The triangles cross where those two stretches
%   overlap, so the answer is the overlap of two intervals on that line.

DoesCross = false;
Start = zeros(1, 3);
Finish = zeros(1, 3);

NormalA = cross(TriangleA(2, :) - TriangleA(1, :), TriangleA(3, :) - TriangleA(1, :));
NormalB = cross(TriangleB(2, :) - TriangleB(1, :), TriangleB(3, :) - TriangleB(1, :));
if norm(NormalA) == 0 || norm(NormalB) == 0
    return   % a triangle with no area cannot cross anything
end
NormalA = NormalA / norm(NormalA);
NormalB = NormalB / norm(NormalB);

% Signed distance of each corner of one triangle from the other's plane.
% Corners within Tolerance of the plane count as exactly on it.
DistanceOfB = (TriangleB - TriangleA(1, :)) * NormalA';
DistanceOfA = (TriangleA - TriangleB(1, :)) * NormalB';
DistanceOfB(abs(DistanceOfB) <= Tolerance) = 0;
DistanceOfA(abs(DistanceOfA) <= Tolerance) = 0;

% Entirely on one side of the other's plane: no crossing.
if all(DistanceOfB > 0) || all(DistanceOfB < 0) || ...
   all(DistanceOfA > 0) || all(DistanceOfA < 0)
    return
end

% Entirely in the other's plane: coplanar, which is not handled.
if all(DistanceOfB == 0) || all(DistanceOfA == 0)
    return
end

% The planes cross along a line with this direction.
LineDirection = cross(NormalA, NormalB);
if norm(LineDirection) < 1e-12
    return   % parallel planes that are not coplanar
end
LineDirection = LineDirection / norm(LineDirection);

% Where each triangle meets the line, as points and positions along it.
[PointsA, PositionsA] = FindPlaneCrossings(TriangleA, DistanceOfA, LineDirection);
[PointsB, PositionsB] = FindPlaneCrossings(TriangleB, DistanceOfB, LineDirection);
if numel(PositionsA) < 2 || numel(PositionsB) < 2
    return   % one of them only touches the plane at a single corner
end

% The overlap runs from the later of the two starts to the earlier of the
% two ends. Each end is a point already known from one of the triangles.
[MinimumA, MinimumIndexA] = min(PositionsA);
[MaximumA, MaximumIndexA] = max(PositionsA);
[MinimumB, MinimumIndexB] = min(PositionsB);
[MaximumB, MaximumIndexB] = max(PositionsB);

if MinimumA >= MinimumB
    OverlapStart = MinimumA;
    Start = PointsA(MinimumIndexA, :);
else
    OverlapStart = MinimumB;
    Start = PointsB(MinimumIndexB, :);
end

if MaximumA <= MaximumB
    OverlapEnd = MaximumA;
    Finish = PointsA(MaximumIndexA, :);
else
    OverlapEnd = MaximumB;
    Finish = PointsB(MaximumIndexB, :);
end

% Overlaps that are a single point (or shorter than Tolerance) do not cut
% anything.
DoesCross = (OverlapEnd - OverlapStart) > Tolerance;

end


function [Points, Positions] = FindPlaneCrossings(Triangle, Distances, LineDirection)
%FINDPLANECROSSINGS Points where a triangle meets a plane, with positions along a line.
%   Distances holds the signed distance of each corner from the plane (0 for
%   corners on it). A corner on the plane is a crossing point itself; an
%   edge whose ends are on opposite sides crosses in between. Positions is
%   the distance of each point along LineDirection, which lets points from
%   different triangles be compared.

Points = zeros(0, 3);

for Corner = 1:3
    NextCorner = mod(Corner, 3) + 1;

    if Distances(Corner) == 0
        Points(end + 1, :) = Triangle(Corner, :); %#ok<AGROW>
    elseif Distances(Corner) * Distances(NextCorner) < 0
        % Opposite sides: the distance reaches zero after the fraction
        % Da / (Da - Db) of the way from this corner to the next.
        Fraction = Distances(Corner) / (Distances(Corner) - Distances(NextCorner));
        Points(end + 1, :) = Triangle(Corner, :) + ...
            Fraction * (Triangle(NextCorner, :) - Triangle(Corner, :)); %#ok<AGROW>
    end
end

Positions = Points * LineDirection';

end


function NewFaces = SplitIntersectedTriangles(Faces, Verts, CutFace, CutSegments, Tolerance)
%SPLITINTERSECTEDTRIANGLES Retriangulate every triangle that has segments in it.
%   CutFace is a Kx1 list of triangle numbers and CutSegments the matching
%   Kx2 end-point vertex indices; a triangle can appear several times.
%   Triangles without segments are returned unchanged.

if isempty(CutFace)
    NewFaces = Faces;
    return
end

% Gather the segments belonging to each cut triangle.
[CutFaceIds, ~, GroupOfSegment] = unique(CutFace);
SegmentsPerFace = accumarray(GroupOfSegment, (1:numel(CutFace))', [], ...
    @(Rows) {CutSegments(Rows, :)});

UncutMask = true(size(Faces, 1), 1);
UncutMask(CutFaceIds) = false;

NewPieces = cell(numel(CutFaceIds), 1);
for Index = 1:numel(CutFaceIds)
    NewPieces{Index} = RetriangulateFace( ...
        Faces(CutFaceIds(Index), :), SegmentsPerFace{Index}, Verts, Tolerance);
end

NewFaces = [Faces(UncutMask, :); vertcat(NewPieces{:})];

end


function Pieces = RetriangulateFace(CornerIds, Segments, Verts, Tolerance)
%RETRIANGULATEFACE Split one triangle so that every segment becomes an edge.
%   CornerIds is the 1x3 triangle, Segments an Kx2 list of vertex index
%   pairs lying inside it. Returns the replacement triangles with the
%   winding of the original.

% Segments that start and end at the same vertex, and repeats of the same
% segment (from several neighboring triangles), add nothing.
Segments = Segments(Segments(:, 1) ~= Segments(:, 2), :);
Segments = unique(sort(Segments, 2), "rows");

% A triangle with a repeated corner has no area to subdivide.
if isempty(Segments) || numel(unique(CornerIds)) < 3
    Pieces = CornerIds;
    return
end

% Every vertex used in this triangle: its corners plus the segment ends.
LocalIds = unique([CornerIds(:); Segments(:)]);
LocalPoints = Verts(LocalIds, :);

% Describe the triangle's plane with a right-handed 2D coordinate system
% (AxisU, AxisV) such that AxisU x AxisV = the triangle's normal. A triangle
% that runs counter-clockwise in (U, V) then keeps the original winding.
Corners = Verts(CornerIds, :);
Normal = cross(Corners(2, :) - Corners(1, :), Corners(3, :) - Corners(1, :));
Normal = Normal / norm(Normal);
AxisU = (Corners(2, :) - Corners(1, :)) / norm(Corners(2, :) - Corners(1, :));
AxisV = cross(Normal, AxisU);

Offsets = LocalPoints - Corners(1, :);
Points2D = [Offsets * AxisU', Offsets * AxisV'];

% Points that belong on a triangle edge (where the intersection curve leaves
% the triangle) are moved exactly onto that edge. Otherwise rounding could
% leave one a hair inside, and the edge would not be split there, giving a
% gap in the surface next to the neighboring triangle.
IsCorner = ismember(LocalIds, CornerIds(:));
[~, CornerRows] = ismember(CornerIds(:), LocalIds);
Corners2D = Points2D(CornerRows, :);
for EdgeStart = 1:3
    EdgeEnd = mod(EdgeStart, 3) + 1;
    EdgeVector = Corners2D(EdgeEnd, :) - Corners2D(EdgeStart, :);
    EdgeLength = norm(EdgeVector);
    EdgeDirection = EdgeVector / EdgeLength;

    FromStart = Points2D - Corners2D(EdgeStart, :);
    AlongEdge = FromStart * EdgeDirection';
    % 2D cross product with the edge direction = signed distance to the edge line.
    AcrossEdge = FromStart(:, 1) * EdgeDirection(2) - FromStart(:, 2) * EdgeDirection(1);

    OnEdge = abs(AcrossEdge) <= Tolerance & AlongEdge > 0 & ...
             AlongEdge < EdgeLength & ~IsCorner;
    Points2D(OnEdge, :) = Corners2D(EdgeStart, :) + AlongEdge(OnEdge) * EdgeDirection;
end

% Triangulate with the segments forced to be edges. All points lie inside
% the triangle, so the triangulation covers exactly the triangle.
[~, SegmentsLocal] = ismember(Segments, LocalIds);
try
    Triangulation = delaunayTriangulation(Points2D, SegmentsLocal);
catch Problem
    error("JoinMeshes:RetriangulationFailed", ...
        "Could not retriangulate a cut triangle (%s). Two intersection segments probably cross inside one triangle.", ...
        Problem.message);
end
if size(Triangulation.Points, 1) ~= numel(LocalIds)
    error("JoinMeshes:RetriangulationFailed", ...
        "Two points of a cut triangle coincide. The mesh probably has duplicate vertices or degenerate triangles.");
end
LocalTriangles = Triangulation.ConnectivityList;

% Make every triangle counter-clockwise in (U, V). Twice the signed 2D area
% is the 2D cross product of two edges; negative means clockwise.
CornerA = Points2D(LocalTriangles(:, 1), :);
CornerB = Points2D(LocalTriangles(:, 2), :);
CornerC = Points2D(LocalTriangles(:, 3), :);
TwiceSignedArea = (CornerB(:, 1) - CornerA(:, 1)) .* (CornerC(:, 2) - CornerA(:, 2)) ...
                - (CornerB(:, 2) - CornerA(:, 2)) .* (CornerC(:, 1) - CornerA(:, 1));
IsClockwise = TwiceSignedArea < 0;
LocalTriangles(IsClockwise, :) = LocalTriangles(IsClockwise, [1 3 2]);

% A point that rounding left a hair inside an edge makes a sliver triangle
% along that edge. Dropping it leaves the edge split at the point, which is
% what is wanted.
OriginalTwiceArea = abs((Corners2D(2, 1) - Corners2D(1, 1)) * (Corners2D(3, 2) - Corners2D(1, 2)) ...
                      - (Corners2D(2, 2) - Corners2D(1, 2)) * (Corners2D(3, 1) - Corners2D(1, 1)));
IsSliver = abs(TwiceSignedArea) < 1e-12 * OriginalTwiceArea;
LocalTriangles(IsSliver, :) = [];

% Back to vertex numbers of the whole mesh. reshape keeps a single triangle as 1x3.
Pieces = reshape(LocalIds(LocalTriangles), size(LocalTriangles));

end


function Faces = RemoveRepeatedIndexFaces(Faces)
%REMOVEREPEATEDINDEXFACES Delete triangles that use the same vertex twice.

HasRepeat = Faces(:, 1) == Faces(:, 2) | Faces(:, 2) == Faces(:, 3) | Faces(:, 1) == Faces(:, 3);
Faces(HasRepeat, :) = [];

end


function Faces = RemoveZeroAreaFaces(Faces, Verts, MinimumTwiceArea)
%REMOVEZEROAREAFACES Delete triangles whose area is below a threshold.
%   MinimumTwiceArea is compared with the length of the cross product of two
%   edges, which is twice the triangle area.

EdgeOne = Verts(Faces(:, 2), :) - Verts(Faces(:, 1), :);
EdgeTwo = Verts(Faces(:, 3), :) - Verts(Faces(:, 1), :);
TwiceArea = vecnorm(cross(EdgeOne, EdgeTwo, 2), 2, 2);
Faces(TwiceArea <= MinimumTwiceArea, :) = [];

end


function Faces = RemoveDuplicateFaces(Faces)
%REMOVEDUPLICATEFACES Resolve triangles that occupy the same place twice.
%   Triangles using the same three vertices come from meshes that share a
%   face. If they all have the same winding, one copy is kept. If their
%   windings oppose, the triangle is a wall between two solids and every
%   copy is removed.

SortedFaces = sort(Faces, 2);
[~, ~, Group] = unique(SortedFaces, "rows");
Group = Group(:);

% A triangle is "even" when its index order is a rotation of the sorted
% order, which means it has the same winding as the sorted triangle.
IsEven = all(Faces == SortedFaces, 2) ...
       | all(Faces == SortedFaces(:, [2 3 1]), 2) ...
       | all(Faces == SortedFaces(:, [3 1 2]), 2);

GroupSize = accumarray(Group, 1);
GroupEvenCount = accumarray(Group, double(IsEven));
HasOpposingWinding = GroupEvenCount > 0 & GroupEvenCount < GroupSize;

[~, FirstRows] = unique(Group, "first");
IsFirstOfGroup = false(size(Faces, 1), 1);
IsFirstOfGroup(FirstRows) = true;

Faces = Faces(IsFirstOfGroup & ~HasOpposingWinding(Group), :);

end


function [SideRegion, FaceComponent] = FindSideRegions(Faces, Verts)
%FINDSIDEREGIONS Label the regions of space the triangles divide, per triangle side.
%   Triangle k has a front side (node 2k-1, where its normal points) and a
%   back side (node 2k). Two sides belong to the same region when they face
%   the same wedge between neighboring triangles around some edge.
%
%   Around an edge, every triangle meeting there extends away from the edge
%   in some direction. Sorting those directions by angle about the edge
%   gives the wedges: the wedge counter-clockwise of one triangle is also
%   the wedge clockwise of the next triangle in the sorted order. Which
%   side of a triangle faces its counter-clockwise wedge depends on how the
%   triangle runs along the edge: if its winding follows the edge from the
%   lower to the higher vertex number, the front side does (by the right-hand
%   rule), otherwise the back side does.
%
%   SideRegion is a 2N x 1 region label per side. FaceComponent is an N x 1
%   label of the group of edge-connected triangles each triangle is in.

FaceCount = size(Faces, 1);
NodeCount = 2 * FaceCount;

% Three half-edges per triangle (corner 1 to 2, 2 to 3, 3 to 1), listed as
% three blocks of FaceCount rows. The third corner is the one off the edge.
HalfStart = [Faces(:, 1); Faces(:, 2); Faces(:, 3)];
HalfEnd = [Faces(:, 2); Faces(:, 3); Faces(:, 1)];
HalfThird = [Faces(:, 3); Faces(:, 1); Faces(:, 2)];
HalfFace = repmat((1:FaceCount)', 3, 1);
HalfCount = 3 * FaceCount;

% Edges are shared between triangles, so each is identified by its two
% vertex numbers in increasing order.
IsForward = HalfStart < HalfEnd;
LowVert = min(HalfStart, HalfEnd);
HighVert = max(HalfStart, HalfEnd);
[~, ~, EdgeId] = unique([LowVert, HighVert], "rows");
EdgeId = EdgeId(:);

% Direction along the edge, and the direction from the edge into the
% triangle (perpendicular to the edge).
Axis = Verts(HighVert, :) - Verts(LowVert, :);
Axis = Axis ./ vecnorm(Axis, 2, 2);
IntoFace = Verts(HalfThird, :) - Verts(LowVert, :);
IntoFace = IntoFace - dot(IntoFace, Axis, 2) .* Axis;

% Angle of each triangle about its edge, measured from the first triangle
% listed for that edge. Angles increase counter-clockwise about Axis.
[~, FirstHalfOfEdge] = unique(EdgeId, "first");
Reference = IntoFace(FirstHalfOfEdge(EdgeId), :);
Reference = Reference ./ vecnorm(Reference, 2, 2);
Perpendicular = cross(Axis, Reference, 2);
Angle = atan2(dot(IntoFace, Perpendicular, 2), dot(IntoFace, Reference, 2));

% Sort by edge, then by angle. The next triangle around the edge is the
% following row, wrapping back to the first row of the edge after the last.
[~, Order] = sortrows([EdgeId, Angle]);
SortedEdge = EdgeId(Order);
SortedFace = HalfFace(Order);
SortedForward = IsForward(Order);

IsFirstOfEdge = [true; SortedEdge(1:end - 1) ~= SortedEdge(2:end)];
IsLastOfEdge = [IsFirstOfEdge(2:end); true];
EdgeNumber = cumsum(IsFirstOfEdge);
FirstRowOfEdge = find(IsFirstOfEdge);

NextRow = (1:HalfCount)' + 1;
NextRow(IsLastOfEdge) = FirstRowOfEdge(EdgeNumber(IsLastOfEdge));

% The side of each triangle facing its counter-clockwise wedge, and the
% side facing its clockwise wedge.
FrontNode = 2 * SortedFace - 1;
BackNode = 2 * SortedFace;
CounterClockwiseNode = BackNode;
CounterClockwiseNode(SortedForward) = FrontNode(SortedForward);
ClockwiseNode = FrontNode;
ClockwiseNode(SortedForward) = BackNode(SortedForward);

% The wedge between a triangle and the next one around the edge touches the
% first one's counter-clockwise side and the next one's clockwise side.
WedgeLinks = sparse(CounterClockwiseNode, ClockwiseNode(NextRow), 1, NodeCount, NodeCount);
SideRegion = conncomp(graph(spones(WedgeLinks + WedgeLinks')));
SideRegion = SideRegion(:);

% Triangles are in the same component when their sides are linked through
% edges or when they are the front and back of one triangle.
FaceIds = (1:FaceCount)';
FaceLinks = sparse(2 * FaceIds - 1, 2 * FaceIds, 1, NodeCount, NodeCount);
AllLinks = WedgeLinks + FaceLinks;
NodeComponent = conncomp(graph(spones(AllLinks + AllLinks')));
FaceComponent = NodeComponent(2 * FaceIds);
FaceComponent = FaceComponent(:);

end


function ExteriorSide = FindExteriorSide(Faces, Verts, FaceIds)
%FINDEXTERIORSIDE A triangle side that faces the outside of a group of triangles.
%   Faces are the triangles of one connected group and FaceIds their numbers
%   in the full list. A ray started far outside the group, traveling toward
%   -x, reaches the group from the outside, so the side of the first
%   triangle it hits that faces the ray start is an outside side. Returns
%   the side's node number (2k-1 for the front of triangle k, 2k for the back).

Corner0 = Verts(Faces(:, 1), :);
Corner1 = Verts(Faces(:, 2), :);
Corner2 = Verts(Faces(:, 3), :);
Edge1 = Corner1 - Corner0;
Edge2 = Corner2 - Corner0;
Normals = cross(Edge1, Edge2, 2);
Centers = (Corner0 + Corner1 + Corner2) / 3;

RayDirection = [-1 0 0];
RayDirections = repmat(RayDirection, size(Faces, 1), 1);
StartX = max(Verts(unique(Faces(:)), 1)) + 1;

% Ray-triangle test (Moller-Trumbore) set up once for all triangles. A
% triangle the ray runs parallel to cannot be hit.
CrossDirectionEdge2 = cross(RayDirections, Edge2, 2);
Determinant = dot(Edge1, CrossDirectionEdge2, 2);
IsParallel = abs(Determinant) <= 1e-12 * vecnorm(Edge1, 2, 2) .* vecnorm(Edge2, 2, 2);

% Aim the ray at the center of one of the triangles farthest in +x, so it is
% sure to reach the group. If the nearest hit is too close to a triangle
% edge to trust, try the next candidate.
[~, ByDescendingX] = sort(Centers(:, 1), "descend");
for Target = ByDescendingX(1:min(50, numel(ByDescendingX)))'
    Origin = [StartX, Centers(Target, 2), Centers(Target, 3)];

    FromCorner = Origin - Corner0;
    U = dot(FromCorner, CrossDirectionEdge2, 2) ./ Determinant;
    CrossFromEdge1 = cross(FromCorner, Edge1, 2);
    V = dot(RayDirections, CrossFromEdge1, 2) ./ Determinant;
    T = dot(Edge2, CrossFromEdge1, 2) ./ Determinant;

    IsHit = ~IsParallel & U >= 0 & V >= 0 & (U + V) <= 1 & T > 0;
    if ~any(IsHit)
        continue
    end

    HitRows = find(IsHit);
    [~, Nearest] = min(T(HitRows));
    HitRow = HitRows(Nearest);

    % Reject hits that graze a triangle edge, where the nearest triangle
    % is not clearly defined.
    if min([U(HitRow), V(HitRow), 1 - U(HitRow) - V(HitRow)]) < 1e-6
        continue
    end

    % The ray travels along RayDirection. If the triangle's normal points
    % against it, the ray meets the front side; otherwise the back side.
    MeetsFront = dot(Normals(HitRow, :), RayDirection) < 0;
    ExteriorSide = 2 * FaceIds(HitRow) - MeetsFront;
    return
end

error("JoinMeshes:ExteriorNotFound", ...
    "Could not tell which side of the surface is the outside. The mesh may be degenerate.");

end


function [Faces, Verts] = RemoveUnusedVerts(Faces, Verts)
%REMOVEUNUSEDVERTS Drop vertices no face uses and renumber the faces.

UsedVerts = unique(Faces(:));
NewIndex = zeros(size(Verts, 1), 1);
NewIndex(UsedVerts) = 1:numel(UsedVerts);

Verts = Verts(UsedVerts, :);
Faces = reshape(NewIndex(Faces), size(Faces));

end
