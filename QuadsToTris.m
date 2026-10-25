function [TriFaces, Verts] = QuadsToTris(Faces, Verts)
%QUADSTOTRIS Convert the quad faces of a mesh into triangles.
%
%   [TriFaces, Verts] = QuadsToTris(Faces, Verts)
%
%   Each quad [A B C D] is split along the A-C diagonal into the two
%   triangles [A B C] and [A C D]. Face order is preserved: the two
%   triangles from a quad appear next to each other where the quad used to
%   be. Meshes are assumed to be all quads or all triangles (not mixed); an
%   all-triangle mesh is returned unchanged.
%
%   Inputs
%       Faces - Nx4 (quads) or Nx3 (triangles) matrix of 1-based vertex
%               indices.
%       Verts - Mx3 matrix of vertex positions.
%
%   Outputs
%       TriFaces - Kx3 matrix of triangle vertex indices (K = 2N for quads).
%       Verts    - Unchanged copy of the input vertices. No vertices are
%                  added or moved, so this is returned only so the pair can
%                  be passed along together.

arguments
    Faces (:,:) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
end

% Faces must have 3 or 4 columns (triangles or quads).
FaceWidth = size(Faces, 2);
if FaceWidth ~= 3 && FaceWidth ~= 4
    error("QuadsToTris:InvalidFaceWidth", ...
        "Faces must have 3 or 4 columns, but has %d.", FaceWidth);
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("QuadsToTris:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

% Nothing to split if the input is already all triangles.
if FaceWidth == 3
    TriFaces = Faces;
    return
end

% Each quad becomes two triangles, so the output has twice as many rows.
QuadCount = size(Faces, 1);
TriFaces = zeros(2 * QuadCount, 3, "like", Faces);

% Odd rows get the first triangle [A B C] of every quad; even rows get the
% second triangle [A C D]. Interleaving keeps the original face order.
TriFaces(1:2:end, :) = Faces(:, [1 2 3]);
TriFaces(2:2:end, :) = Faces(:, [1 3 4]);

end
