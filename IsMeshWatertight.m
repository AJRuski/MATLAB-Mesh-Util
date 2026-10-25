function [IsWatertight, OpenEdges, NonManifoldEdges] = IsMeshWatertight(Faces, Verts, Options)
%ISMESHWATERTIGHT Check whether a triangle mesh is watertight.
%
%   IsWatertight = IsMeshWatertight(Faces, Verts)
%   [IsWatertight, OpenEdges, NonManifoldEdges] = IsMeshWatertight(Faces, Verts)
%   ___ = IsMeshWatertight(Faces, Verts, ShowFigure=true)
%
%   A mesh is watertight when every edge is shared by exactly two
%   triangles. Edges are compared by their vertex indices only, so
%   duplicated vertices along a seam (same position, different index) show
%   up as open edges.
%
%   Not checked: face winding direction (see FixFaceNormals) and vertices
%   where two surfaces touch at a single point.
%
%   Inputs
%       Faces - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts - Mx3 matrix of vertex positions.
%
%   Name-value options
%       ShowFigure - If true and the mesh is not watertight, draw the mesh
%                    with the problem edges highlighted. Default false.
%
%   Outputs
%       IsWatertight     - True if every edge is used by exactly two triangles.
%       OpenEdges        - Kx2 vertex index pairs for edges used by only one
%                          triangle (holes or gaps).
%       NonManifoldEdges - Jx2 vertex index pairs for edges used by more than
%                          two triangles.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
    Options.ShowFigure (1,1) logical = false
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("IsMeshWatertight:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

% An empty mesh is not watertight, and there are no edges to report.
if isempty(Faces)
    IsWatertight = false;
    OpenEdges = zeros(0, 2);
    NonManifoldEdges = zeros(0, 2);
    return
end

% List every edge of every triangle. A triangle [A B C] has the edges
% A->B, B->C and C->A.
EdgeStart = Faces(:);
EdgeEnd = reshape(Faces(:, [2 3 1]), [], 1);

% Sort each edge's endpoints so A->B and B->A count as the same edge, then
% find the unique edges and how many triangles use each one.
[UniqueEdges, ~, EdgeId] = unique(sort([EdgeStart, EdgeEnd], 2), "rows");
UsesPerEdge = accumarray(EdgeId, 1);

% Edges used once are boundaries of a hole; edges used more than twice are
% non-manifold (three or more triangles meeting along one edge).
OpenEdges = UniqueEdges(UsesPerEdge == 1, :);
NonManifoldEdges = UniqueEdges(UsesPerEdge > 2, :);

% Watertight means no problem edges at all.
IsWatertight = isempty(OpenEdges) && isempty(NonManifoldEdges);

% Only draw when asked, and only when there is something to show.
if Options.ShowFigure && ~IsWatertight
    ShowProblemEdges(Faces, Verts, OpenEdges, NonManifoldEdges);
end

end


function ShowProblemEdges(Faces, Verts, OpenEdges, NonManifoldEdges)
%SHOWPROBLEMEDGES Draw the mesh with open and non-manifold edges highlighted.

figure("Name", "Mesh watertight check");
hold on

% Draw the mesh as a translucent surface so the highlighted edges stand out.
patch("Faces", Faces, "Vertices", Verts, ...
    "FaceColor", [0.7 0.7 0.7], "FaceAlpha", 0.4, "EdgeColor", "none");

% Highlight each kind of problem edge, and build the legend only from the
% kinds that actually exist.
LegendHandles = gobjects(0);
LegendNames = strings(0);

if ~isempty(OpenEdges)
    LegendHandles(end+1) = PlotEdges(Verts, OpenEdges, "r");
    LegendNames(end+1) = sprintf("Open edges (%d)", size(OpenEdges, 1));
end
if ~isempty(NonManifoldEdges)
    LegendHandles(end+1) = PlotEdges(Verts, NonManifoldEdges, "m");
    LegendNames(end+1) = sprintf("Non-manifold edges (%d)", size(NonManifoldEdges, 1));
end

legend(LegendHandles, LegendNames);
title("Mesh is not watertight");
axis equal
view(3)
grid on
hold off

end


function LineHandle = PlotEdges(Verts, Edges, Color)
%PLOTEDGES Plot many edges as one line object.
%   Each edge becomes "start point, end point, NaN". The NaN breaks the line
%   so all edges can be drawn with a single plot3 call.

Starts = Verts(Edges(:, 1), :);
Ends = Verts(Edges(:, 2), :);
Gaps = NaN(size(Edges, 1), 1);

% Interleave start, end, gap per edge, then flatten each coordinate into a
% single column.
X = [Starts(:, 1), Ends(:, 1), Gaps]';
Y = [Starts(:, 2), Ends(:, 2), Gaps]';
Z = [Starts(:, 3), Ends(:, 3), Gaps]';

LineHandle = plot3(X(:), Y(:), Z(:), Color, "LineWidth", 2.5);

end
