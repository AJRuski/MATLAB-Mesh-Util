function [FixedFaces, FixedVerts, RemovedFaceCount] = FixDegenerateTriangles(Faces, Verts, Options)
%FIXDEGENERATETRIANGLES Remove degenerate triangles by collapsing their shortest edge.
%
%   [FixedFaces, FixedVerts] = FixDegenerateTriangles(Faces, Verts)
%   [FixedFaces, FixedVerts, RemovedFaceCount] = FixDegenerateTriangles(Faces, Verts)
%   ___ = FixDegenerateTriangles(Faces, Verts, Tolerance=1e-8, MaxPasses=50)
%
%   Degenerate triangles are found with FindDegenerateTriangles and fixed
%   in two ways:
%     - Triangles that repeat a vertex index (e.g. [4 4 9]) are simply
%       deleted. For a watertight mesh this is safe: such a triangle adds
%       an extra pair of uses to one edge, so removing it restores the
%       edge to its normal two uses.
%     - Triangles with zero area but three different vertices are fixed by
%       an edge collapse: the shortest edge is shrunk to a point by merging
%       its two vertices at the midpoint of the edge. The degenerate
%       triangle, and the triangle on the other side of that edge, then
%       have a repeated vertex and are deleted. The neighbors are rejoined
%       along the merged vertex, so a watertight mesh stays watertight in
%       the usual case.
%
%   Collapsing one edge can create new degenerate triangles nearby, so the
%   process repeats until none remain or MaxPasses is reached.
%
%   Vertices that no face uses are dropped at the end and the face indices
%   are renumbered, so FixedVerts can have fewer rows than Verts and its
%   row numbers will not match the original vertex numbers.
%
%   The result is not guaranteed to stay watertight (an edge collapse can
%   occasionally pinch the surface). If the input was watertight and the
%   output is not, a warning is issued. Run FixFaceNormals afterwards if you
%   need to confirm the normals are still consistent.
%
%   Inputs
%       Faces - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts - Mx3 matrix of vertex positions.
%
%   Name-value options
%       Tolerance - Relative area tolerance passed to
%                   FindDegenerateTriangles. Default 1e-10.
%       MaxPasses - Maximum number of find-and-fix rounds. Default 50.
%
%   Outputs
%       FixedFaces       - Kx3 triangle vertex indices with degenerate
%                          triangles removed.
%       FixedVerts       - Px3 vertex positions, with unused vertices removed.
%       RemovedFaceCount - Number of triangles deleted (N - K).

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
    Options.Tolerance (1,1) double {mustBeNonnegative, mustBeFinite} = 1e-10
    Options.MaxPasses (1,1) double {mustBeInteger, mustBePositive} = 50
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("FixDegenerateTriangles:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

% Remember the starting state so we can report the change and check whether
% watertightness was lost.
OriginalFaceCount = size(Faces, 1);
WasWatertight = IsMeshWatertight(Faces, Verts);

% Work in double precision since vertices are moved to edge midpoints.
Verts = double(Verts);

%% Fix in rounds until no degenerate triangles are left

for Pass = 1:Options.MaxPasses

    % Delete triangles that reuse a vertex index. Both the originals and
    % any created by collapses in the previous round end up here.
    [~, HasRepeatedIndex] = FindDegenerateTriangles(Faces, Verts, ...
        Tolerance=Options.Tolerance);
    Faces(HasRepeatedIndex, :) = [];

    % Find the triangles that still have (nearly) zero area.
    IsDegenerate = FindDegenerateTriangles(Faces, Verts, ...
        Tolerance=Options.Tolerance);
    if ~any(IsDegenerate)
        break   % mesh is clean
    end

    % Collapse the shortest edge of each degenerate triangle, one at a time.
    % Earlier collapses in this loop may already have changed a later
    % triangle, so its current corners are looked up fresh each time.
    for FaceIndex = find(IsDegenerate)'
        Corners = Faces(FaceIndex, :);

        % If an earlier collapse already merged two of this triangle's
        % corners, it is handled by the deletion step at the top of the
        % next round.
        if numel(unique(Corners)) < 3
            continue
        end

        % Pick the shortest of the three edges A-B, B-C and C-A.
        EdgeStarts = Corners([1 2 3]);
        EdgeEnds = Corners([2 3 1]);
        EdgeLengths = vecnorm(Verts(EdgeEnds, :) - Verts(EdgeStarts, :), 2, 2);
        [~, ShortestEdge] = min(EdgeLengths);
        KeptVertex = EdgeStarts(ShortestEdge);
        RemovedVertex = EdgeEnds(ShortestEdge);

        % Merge the two vertices: the kept vertex moves to the midpoint,
        % and every face that used the removed vertex now uses the kept one.
        Verts(KeptVertex, :) = (Verts(KeptVertex, :) + Verts(RemovedVertex, :)) / 2;
        Faces(Faces == RemovedVertex) = KeptVertex;
    end
end

% Final cleanup of the triangles that collapsed in the last round.
[~, HasRepeatedIndex] = FindDegenerateTriangles(Faces, Verts, ...
    Tolerance=Options.Tolerance);
Faces(HasRepeatedIndex, :) = [];

% If rounds ran out before the mesh was clean, say so.
if any(FindDegenerateTriangles(Faces, Verts, Tolerance=Options.Tolerance))
    warning("FixDegenerateTriangles:NotFullyFixed", ...
        "Degenerate triangles remain after %d passes.", Options.MaxPasses);
end

%% Drop unused vertices and renumber the faces

% UsedVerts lists the vertices still referenced by some face. NewIndex maps
% each old vertex number to its new position in that list.
UsedVerts = unique(Faces(:));
NewIndex = zeros(size(Verts, 1), 1);
NewIndex(UsedVerts) = 1:numel(UsedVerts);

FixedVerts = Verts(UsedVerts, :);
FixedFaces = reshape(NewIndex(Faces), size(Faces));   % reshape keeps a single face as 1x3

RemovedFaceCount = OriginalFaceCount - size(FixedFaces, 1);

% An edge collapse can pinch a closed surface; warn if that happened.
if WasWatertight && ~IsMeshWatertight(FixedFaces, FixedVerts)
    warning("FixDegenerateTriangles:LostWatertight", ...
        "The input was watertight but the result is not. Inspect it with IsMeshWatertight(..., ShowFigure=true).");
end

end
