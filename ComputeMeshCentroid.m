function [RelativeCentroid, Centroid] = ComputeMeshCentroid(Faces, Verts, ReferencePoint)
%COMPUTEMESHCENTROID Compute the volume centroid of a closed mesh relative to a point.
%
%   RelativeCentroid = ComputeMeshCentroid(Faces, Verts, ReferencePoint)
%   [RelativeCentroid, Centroid] = ComputeMeshCentroid(Faces, Verts, ReferencePoint)
%
%   The centroid is the center of mass of the solid enclosed by the mesh,
%   assuming uniform density. It is found by treating the mesh as a set of
%   tetrahedra, each formed by one triangle and a common origin point. Each
%   tetrahedron has a signed volume of dot(A, cross(B, C)) / 6 and its own
%   centroid at (A + B + C) / 4 (the origin corner contributes zero). The
%   mesh centroid is the volume-weighted average of those centroids.
%
%   The mesh is assumed to be watertight. No checks are made for this.
%   Because the volume appears in both the numerator and denominator, the
%   result is the same whether the face normals point outward or inward.
%
%   Inputs
%       Faces          - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts          - Mx3 matrix of vertex positions.
%       ReferencePoint - 1x3 [x y z] point to measure from. Default [0 0 0].
%
%   Outputs
%       RelativeCentroid - 1x3 vector from ReferencePoint to the centroid
%                          (Centroid - ReferencePoint).
%       Centroid         - 1x3 centroid in the same coordinates as Verts.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
    ReferencePoint (1,3) {mustBeNumeric, mustBeReal, mustBeFinite} = [0 0 0]
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("ComputeMeshCentroid:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

% Work in double precision, with the mesh shifted so the average vertex is
% at the origin. The centroid does not depend on the origin used for the
% tetrahedra, but a mesh far from the origin would otherwise lose
% precision. The shift is added back at the end.
Verts = double(Verts);
ShiftToOrigin = mean(Verts, 1);
Verts = Verts - ShiftToOrigin;

% Gather the three corners of every triangle as Nx3 arrays.
CornerA = Verts(Faces(:, 1), :);
CornerB = Verts(Faces(:, 2), :);
CornerC = Verts(Faces(:, 3), :);

% Signed volume of the tetrahedron (origin, A, B, C) for each triangle.
TetVolumes = dot(CornerA, cross(CornerB, CornerC, 2), 2) / 6;

% Centroid of each tetrahedron: the average of its four corners, one of
% which is the origin (0), so the sum is just A + B + C divided by 4.
TetCentroids = (CornerA + CornerB + CornerC) / 4;

% The mesh needs a nonzero enclosed volume for the average to make sense.
TotalVolume = sum(TetVolumes);
if TotalVolume == 0
    error("ComputeMeshCentroid:ZeroVolume", ...
        "The mesh encloses zero volume, so its centroid is undefined.");
end

% Volume-weighted average of the tetrahedron centroids (1x3), then undo
% the shift so the result is in the original coordinates.
Centroid = sum(TetVolumes .* TetCentroids, 1) / TotalVolume + ShiftToOrigin;

% Vector from the reference point to the centroid.
RelativeCentroid = Centroid - ReferencePoint;

end
