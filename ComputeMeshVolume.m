function Volume = ComputeMeshVolume(Faces, Verts)
%COMPUTEMESHVOLUME Compute the enclosed volume of a closed mesh.
%
%   Volume = ComputeMeshVolume(Faces, Verts)
%
%   Uses the divergence theorem: the mesh is treated as a set of tetrahedra,
%   each formed by one triangle and the origin. The signed volume of each
%   tetrahedron is dot(A, cross(B, C)) / 6, and the sum over all triangles
%   is the volume enclosed by the surface.
%
%   The mesh is assumed to be watertight with consistently outward-facing
%   face normals (counter-clockwise winding seen from outside). No checks
%   are made for this. If the normals point inward, the result is negative.
%
%   Inputs
%       Faces - Nx3 matrix of 1-based vertex indices (triangles). Use
%               QuadsToTris first if the mesh has quads.
%       Verts - Mx3 matrix of vertex positions.
%
%   Outputs
%       Volume - Scalar signed volume, in the cubed units of Verts.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("ComputeMeshVolume:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

% Work in double precision, and shift the mesh so its centroid sits at the
% origin. The volume does not depend on the origin, but a mesh far from the
% origin would otherwise sum large, nearly cancelling terms and lose
% precision.
Verts = double(Verts);
Verts = Verts - mean(Verts, 1);

% Gather the three corners of every triangle as Nx3 arrays.
CornerA = Verts(Faces(:, 1), :);
CornerB = Verts(Faces(:, 2), :);
CornerC = Verts(Faces(:, 3), :);

% Signed volume of the tetrahedron (origin, A, B, C) for each triangle:
% the scalar triple product A . (B x C), divided by 6.
TetVolumes = dot(CornerA, cross(CornerB, CornerC, 2), 2) / 6;

% The mesh volume is the sum of all the tetrahedron volumes.
Volume = sum(TetVolumes);

end
