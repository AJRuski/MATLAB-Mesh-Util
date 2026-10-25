function [FixedFaces, Verts, FlipMask] = FixFaceNormals(Faces, Verts)
%FIXFACENORMALS Make all triangle normals consistent and pointing outward.
%
%   [FixedFaces, Verts, FlipMask] = FixFaceNormals(Faces, Verts)
%
%   Works in two steps:
%     1. Consistency. Two neighboring triangles are consistent when they
%        travel along their shared edge in opposite directions. Starting
%        from one triangle of each connected piece of the mesh, the
%        orientation is spread across shared edges, flipping any triangle
%        that disagrees with the triangle it was reached from.
%     2. Outward. The signed volume of each connected piece is computed. A
%        negative volume means the piece is inside-out, so all of its
%        triangles are flipped.
%
%   A triangle is flipped by swapping its last two vertices
%   ([A B C] -> [A C B]).
%
%   The mesh must be watertight (every edge shared by exactly two
%   triangles); an error is raised otherwise. Use QuadsToTris first if the
%   mesh has quads.
%
%   Limitation: each connected piece is oriented to point away from its own
%   enclosed volume. A separate inner shell that bounds a hollow cavity
%   would therefore point into the material instead of into the cavity.
%
%   Inputs
%       Faces - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts - Mx3 matrix of vertex positions.
%
%   Outputs
%       FixedFaces - Faces with consistent, outward-facing winding.
%       Verts      - Unchanged copy of the input vertices.
%       FlipMask   - Nx1 logical, true for each triangle that was flipped.
%                    any(FlipMask) is false if the input was already
%                    correct.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("FixFaceNormals:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

FaceCount = size(Faces, 1);

% An empty mesh has nothing to fix.
if FaceCount == 0
    FixedFaces = Faces;
    FlipMask = false(0, 1);
    return
end

%% Step 1: find every edge and which two faces share it

% List every directed edge of every triangle. Row k of EdgeStart/EdgeEnd is
% one edge, and FaceOfEdge says which triangle it belongs to. A triangle
% [A B C] has the edges A->B, B->C and C->A.
EdgeStart = Faces(:);
EdgeEnd = reshape(Faces(:, [2 3 1]), [], 1);
FaceOfEdge = repmat((1:FaceCount)', 3, 1);

% Give the two directions of an edge the same ID by sorting its endpoints,
% then count how many triangle-edges use each ID.
[~, ~, EdgeId] = unique(sort([EdgeStart, EdgeEnd], 2), "rows");
UsesPerEdge = accumarray(EdgeId, 1);

% The mesh is assumed watertight, so every edge must be used exactly twice.
% Stop with a clear message if not, since the pairing below relies on it.
if any(UsesPerEdge ~= 2)
    error("FixFaceNormals:NotWatertight", ...
        "Mesh is not watertight: %d edges are not shared by exactly two faces.", ...
        nnz(UsesPerEdge ~= 2));
end

% Put the two uses of each edge next to each other (sorted by edge ID),
% then split them into a "first" and "second" use.
[~, SortOrder] = sort(EdgeId);
FirstUse = SortOrder(1:2:end);
SecondUse = SortOrder(2:2:end);

% The two triangles touching each edge, and whether they disagree.
% They disagree (Mismatch = true) when both travel the edge in the same
% direction, since consistent neighbors travel it in opposite directions.
FaceA = FaceOfEdge(FirstUse);
FaceB = FaceOfEdge(SecondUse);
Mismatch = EdgeStart(FirstUse) == EdgeStart(SecondUse);

%% Step 2: build a neighbor list for every face

% Store each neighbor relation in both directions, sorted by source face,
% so that the neighbors of face F are rows NeighborStart(F) to
% NeighborStart(F + 1) - 1.
Source = [FaceA; FaceB];
Neighbor = [FaceB; FaceA];
NeighborMismatch = [Mismatch; Mismatch];
[Source, SortOrder] = sort(Source);
Neighbor = Neighbor(SortOrder);
NeighborMismatch = NeighborMismatch(SortOrder);
NeighborStart = [1; 1 + cumsum(accumarray(Source, 1, [FaceCount, 1]))];

%% Step 3: spread the orientation across each connected piece (BFS)

Visited = false(FaceCount, 1);
FlipForConsistency = false(FaceCount, 1);   % flip needed to match the piece's seed face
Component = zeros(FaceCount, 1);            % which connected piece each face is in
ComponentCount = 0;
Queue = zeros(FaceCount, 1);

for Seed = 1:FaceCount
    if Visited(Seed)
        continue   % already reached from an earlier seed
    end

    % Start a new piece. The seed face keeps its orientation (no flip).
    ComponentCount = ComponentCount + 1;
    Visited(Seed) = true;
    Component(Seed) = ComponentCount;
    Queue(1) = Seed;
    Head = 1;
    Tail = 1;

    while Head <= Tail
        Current = Queue(Head);
        Head = Head + 1;

        % Visit each unvisited neighbor of the current face. The neighbor
        % needs a flip relative to the current face when they mismatch, on
        % top of whatever flip the current face already has.
        for Row = NeighborStart(Current):NeighborStart(Current + 1) - 1
            Next = Neighbor(Row);
            if ~Visited(Next)
                Visited(Next) = true;
                Component(Next) = ComponentCount;
                FlipForConsistency(Next) = xor(FlipForConsistency(Current), ...
                    NeighborMismatch(Row));
                Tail = Tail + 1;
                Queue(Tail) = Next;
            end
        end
    end
end

% Sanity check: after the flips, every shared edge must agree. If not, the
% surface cannot be consistently oriented (e.g. a Mobius strip).
StillMismatched = xor(Mismatch, xor(FlipForConsistency(FaceA), FlipForConsistency(FaceB)));
if any(StillMismatched)
    error("FixFaceNormals:NotOrientable", ...
        "The mesh cannot be consistently oriented (%d edges still disagree).", ...
        nnz(StillMismatched));
end

% Apply the consistency flips by swapping the last two vertices.
FixedFaces = Faces;
FixedFaces(FlipForConsistency, :) = FixedFaces(FlipForConsistency, [1 3 2]);

%% Step 4: flip any inside-out piece so it points outward

% Shift to the centroid and use double precision to avoid losing precision
% on meshes far from the origin. The volume does not depend on the origin.
CenteredVerts = double(Verts) - mean(double(Verts), 1);

% Signed volume of the tetrahedron between each triangle and the origin:
% the scalar triple product A . (B x C) / 6.
CornerA = CenteredVerts(FixedFaces(:, 1), :);
CornerB = CenteredVerts(FixedFaces(:, 2), :);
CornerC = CenteredVerts(FixedFaces(:, 3), :);
FaceVolumes = dot(CornerA, cross(CornerB, CornerC, 2), 2) / 6;

% Summing over a closed piece gives its enclosed volume: positive when the
% triangles point outward, negative when the piece is inside-out.
ComponentVolume = accumarray(Component, FaceVolumes, [ComponentCount, 1]);

% Every triangle in a piece with negative volume gets flipped.
FlipForOutward = ComponentVolume(Component) < 0;
FixedFaces(FlipForOutward, :) = FixedFaces(FlipForOutward, [1 3 2]);

% A triangle counts as flipped overall if exactly one of the two steps
% flipped it (two flips cancel out).
FlipMask = xor(FlipForConsistency, FlipForOutward);

end
