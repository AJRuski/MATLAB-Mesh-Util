function [IsDegenerate, HasRepeatedIndex] = FindDegenerateTriangles(Faces, Verts, Options)
%FINDDEGENERATETRIANGLES Find triangles that have no area.
%
%   IsDegenerate = FindDegenerateTriangles(Faces, Verts)
%   [IsDegenerate, HasRepeatedIndex] = FindDegenerateTriangles(Faces, Verts)
%   ___ = FindDegenerateTriangles(Faces, Verts, Tolerance=1e-8)
%   ___ = FindDegenerateTriangles(Faces, Verts, ShowFigure=true)
%
%   A triangle is degenerate when either:
%     - it uses the same vertex index more than once (e.g. [4 4 9]), or
%     - its three corners lie on a line or at a single point, so its area
%       is zero. This covers duplicate vertices at the same position and
%       triangles squashed flat into a line.
%
%   The area test is relative to the triangle's size so it works the same
%   for tiny and huge meshes: a triangle is flagged when
%       2 * Area <= Tolerance * (longest edge length)^2
%   The left side is the length of cross(B - A, C - A); the right side is
%   the area-scale of the triangle's longest edge. The ratio is about the
%   sine of the triangle's smallest angle, so very thin slivers are
%   flagged too.
%
%   Inputs
%       Faces - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts - Mx3 matrix of vertex positions.
%
%   Name-value options
%       Tolerance - Relative area tolerance described above. Use 0 to flag
%                   only exactly zero-area triangles. Default 1e-10.
%       ShowFigure - If true and degenerate triangles exist, draw the mesh
%                    with the degenerate triangles highlighted. Default false.
%
%   Outputs
%       IsDegenerate     - Nx1 logical, true for each degenerate triangle.
%                          Use find(IsDegenerate) for the face numbers.
%       HasRepeatedIndex - Nx1 logical, true for the subset of degenerate
%                          triangles that repeat a vertex index.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
    Options.Tolerance (1,1) double {mustBeNonnegative, mustBeFinite} = 1e-10
    Options.ShowFigure (1,1) logical = false
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("FindDegenerateTriangles:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

%% Triangles that reuse a vertex index

% Compare each pair of the three corners: A-B, B-C and C-A.
HasRepeatedIndex = Faces(:, 1) == Faces(:, 2) ...
                 | Faces(:, 2) == Faces(:, 3) ...
                 | Faces(:, 3) == Faces(:, 1);

%% Triangles with (nearly) zero area

Verts = double(Verts);
CornerA = Verts(Faces(:, 1), :);
CornerB = Verts(Faces(:, 2), :);
CornerC = Verts(Faces(:, 3), :);

% Edge vectors around the triangle.
EdgeAB = CornerB - CornerA;
EdgeBC = CornerC - CornerB;
EdgeCA = CornerA - CornerC;

% Twice the triangle's area is the length of the cross product of two edges.
TwiceArea = vecnorm(cross(EdgeAB, CornerC - CornerA, 2), 2, 2);

% Squared length of the longest edge, used to scale the tolerance to the
% size of the triangle.
LongestEdgeSquared = max([sum(EdgeAB .^ 2, 2), ...
                          sum(EdgeBC .^ 2, 2), ...
                          sum(EdgeCA .^ 2, 2)], [], 2);

% Flag triangles whose area is negligible compared with their size. A
% triangle with all three corners at one point has both sides equal to 0,
% so it is flagged too.
HasZeroArea = TwiceArea <= Options.Tolerance * LongestEdgeSquared;

%% Combine both tests

IsDegenerate = HasRepeatedIndex | HasZeroArea;

% Only draw when asked, and only when there is something to show.
if Options.ShowFigure && any(IsDegenerate)
    ShowDegenerateTriangles(Faces, Verts, IsDegenerate);
end

end


function ShowDegenerateTriangles(Faces, Verts, IsDegenerate)
%SHOWDEGENERATETRIANGLES Draw the mesh with degenerate triangles highlighted.
%   Degenerate triangles have no area, so a filled patch would be
%   invisible. Instead each one is drawn as a red outline plus a red dot at
%   its center (the dot is what shows triangles collapsed to a point).

figure("Name", "Degenerate triangles");
hold on

% Translucent gray mesh for context.
patch("Faces", Faces, "Vertices", Verts, ...
    "FaceColor", [0.7 0.7 0.7], "FaceAlpha", 0.4, "EdgeColor", "none");

% Corners of only the degenerate triangles.
BadFaces = Faces(IsDegenerate, :);
CornerA = Verts(BadFaces(:, 1), :);
CornerB = Verts(BadFaces(:, 2), :);
CornerC = Verts(BadFaces(:, 3), :);
Gaps = NaN(size(BadFaces, 1), 1);

% Outline each triangle as A -> B -> C -> A, followed by a NaN so all
% outlines can be drawn with one plot3 call. Each column of the 5xK arrays
% is one triangle; (:) flattens them into a single line.
X = [CornerA(:, 1), CornerB(:, 1), CornerC(:, 1), CornerA(:, 1), Gaps]';
Y = [CornerA(:, 2), CornerB(:, 2), CornerC(:, 2), CornerA(:, 2), Gaps]';
Z = [CornerA(:, 3), CornerB(:, 3), CornerC(:, 3), CornerA(:, 3), Gaps]';
OutlineHandle = plot3(X(:), Y(:), Z(:), "r-", "LineWidth", 2.5);

% Red dot at the center of each degenerate triangle.
Centers = (CornerA + CornerB + CornerC) / 3;
plot3(Centers(:, 1), Centers(:, 2), Centers(:, 3), "ro", ...
    "MarkerFaceColor", "r", "HandleVisibility", "off");

legend(OutlineHandle, sprintf("Degenerate triangles (%d)", size(BadFaces, 1)));
title("Degenerate triangles");
axis equal
view(3)
grid on
hold off

end
