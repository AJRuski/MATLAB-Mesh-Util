function FigureHandle = VisualizeFaceNormals(Faces, Verts)
%VISUALIZEFACENORMALS Interactively color a mesh by whether its normals face the camera.
%
%   FigureHandle = VisualizeFaceNormals(Faces, Verts)
%
%   Opens a figure showing the mesh. Faces whose normal points toward the
%   camera (the front of the face is what you are looking at) are green;
%   faces whose normal points away from the camera are red. Rotate, pan or
%   zoom the view with the mouse and the colors update as the camera moves.
%
%   For a watertight mesh with correct outward normals, the surface you can
%   see is all green. Red patches that show through indicate inverted
%   normals (or a view into a concavity that curves away from the camera).
%
%   Face normals follow the right-hand rule: for a triangle [A B C] the
%   normal points along cross(B - A, C - A), i.e. counter-clockwise winding
%   as seen from the front.
%
%   Inputs
%       Faces - Nx3 matrix of 1-based vertex indices (triangles).
%       Verts - Mx3 matrix of vertex positions.
%
%   Outputs
%       FigureHandle - Handle to the figure that was created.

arguments
    Faces (:,3) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
    Verts (:,3) {mustBeNumeric, mustBeReal}
end

% Every index must point at an existing vertex.
if any(Faces(:) > size(Verts, 1))
    error("VisualizeFaceNormals:InvalidFaceIndex", ...
        "Faces references a vertex index larger than the %d vertices provided.", ...
        size(Verts, 1));
end

FaceCount = size(Faces, 1);
FacingColor = [0.10 0.70 0.20];   % green: normal points toward the camera
AwayColor = [0.85 0.10 0.10];     % red: normal points away from the camera

%% Precompute each face's center and normal (these never change)

Verts = double(Verts);
CornerA = Verts(Faces(:, 1), :);
CornerB = Verts(Faces(:, 2), :);
CornerC = Verts(Faces(:, 3), :);

FaceCenters = (CornerA + CornerB + CornerC) / 3;

% Only the direction of the normal matters for the facing test, so it is
% left unnormalized.
FaceNormals = cross(CornerB - CornerA, CornerC - CornerA, 2);

%% Build the figure

FigureHandle = figure("Name", "Face normal viewer");
Ax = axes(FigureHandle);

% One patch for the whole mesh with a color per face (flat shading, no
% lighting) so the green/red coloring is not altered by lights.
PatchHandle = patch(Ax, "Faces", Faces, "Vertices", Verts, ...
    "FaceColor", "flat", "FaceVertexCData", zeros(FaceCount, 3), ...
    "EdgeColor", "none");

axis(Ax, "equal");
axis(Ax, "vis3d");   % keep proportions fixed while rotating
view(Ax, 3);
grid(Ax, "on");
title(Ax, "Green = normal faces camera, Red = normal faces away");
rotate3d(FigureHandle, "on");   % left-drag rotates the view

% Color the mesh for the starting camera, then recolor whenever the axes
% finish redrawing. A redraw happens each time the camera moves, which
% gives the real-time update while dragging.
LastIsFacing = [];
UpdateColors();
CameraListener = addlistener(Ax, "MarkedClean", @(~, ~) UpdateColors());

% Store the listener on the figure so it lives as long as the figure does.
setappdata(FigureHandle, "CameraListener", CameraListener);

    function UpdateColors()
        %UPDATECOLORS Recolor faces for the current camera position.

        % Direction from the surface toward the camera. A perspective
        % camera has a different direction for every face; an orthographic
        % camera looks the same way at every face.
        if strcmp(Ax.Projection, "perspective")
            ToCamera = Ax.CameraPosition - FaceCenters;
        else
            ToCamera = Ax.CameraPosition - Ax.CameraTarget;
        end

        % A face points toward the camera when its normal and the
        % direction to the camera make an angle under 90 degrees, i.e.
        % their dot product is positive.
        IsFacing = sum(FaceNormals .* ToCamera, 2) > 0;

        % Changing the colors triggers another redraw (and this callback
        % again), so skip the update when nothing changed.
        if isequal(IsFacing, LastIsFacing)
            return
        end
        LastIsFacing = IsFacing;

        % Pick green for facing faces and red for the rest (Nx1 mask
        % expands against the 1x3 color difference to give an Nx3 array).
        PatchHandle.FaceVertexCData = AwayColor + IsFacing .* (FacingColor - AwayColor);
    end

end
