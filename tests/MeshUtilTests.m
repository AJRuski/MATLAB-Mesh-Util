classdef MeshUtilTests < matlab.unittest.TestCase
%MESHUTILTESTS Unit tests for the mesh utility functions.
%
%   Run from the repository root with:
%       results = runtests("tests");
%       table(results)
%
%   The tests build simple shapes with known answers (box, cylinder,
%   sphere and an L-shaped prism) and compare the library's volume,
%   centroid and plane-cutting results against those answers.
%
%   Known answers used:
%     Box        - exact: volume = width*depth*height, centroid = center.
%     Cylinder   - exact for the N-sided prism that the mesh really is
%                  (area = N/2 * R^2 * sin(2*pi/N)), plus a loose check
%                  against the true cylinder pi*R^2*H.
%     Sphere     - analytic sphere and spherical-cap formulas, with a loose
%                  tolerance because the mesh is a polyhedron approximating
%                  a sphere.
%     L prism    - exact, from the L-shaped outline's area and centroid.

    methods (TestClassSetup)
        function AddRepoToPath(testCase)
            % The tests live in <repo>/tests; the functions are in <repo>.
            RepoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(RepoRoot));
        end
    end

    methods (TestMethodSetup)
        function HideFigures(testCase)
            % Keep figures from popping up during the run, and clean up.
            OriginalVisible = get(groot, "DefaultFigureVisible");
            set(groot, "DefaultFigureVisible", "off");
            testCase.addTeardown(@() set(groot, "DefaultFigureVisible", OriginalVisible));
            testCase.addTeardown(@() close("all", "force"));
        end
    end

    %% QuadsToTris

    methods (Test)
        function QuadsToTrisSingleQuad(testCase)
            Verts = [0 0 0; 1 0 0; 1 1 0; 0 1 0];
            [TriFaces, OutVerts] = QuadsToTris([1 2 3 4], Verts);

            testCase.verifyEqual(TriFaces, [1 2 3; 1 3 4]);
            testCase.verifyEqual(OutVerts, Verts);
        end

        function QuadsToTrisKeepsFaceOrder(testCase)
            Verts = rand(8, 3);
            [TriFaces, ~] = QuadsToTris([1 2 3 4; 5 6 7 8], Verts);

            testCase.verifyEqual(TriFaces, [1 2 3; 1 3 4; 5 6 7; 5 7 8]);
        end

        function QuadsToTrisPreservesArea(testCase)
            % A unit square split into two triangles still has area 1.
            Verts = [0 0 0; 1 0 0; 1 1 0; 0 1 0];
            [TriFaces, ~] = QuadsToTris([1 2 3 4], Verts);

            Normals = cross(Verts(TriFaces(:, 2), :) - Verts(TriFaces(:, 1), :), ...
                            Verts(TriFaces(:, 3), :) - Verts(TriFaces(:, 1), :), 2);
            testCase.verifyEqual(sum(vecnorm(Normals, 2, 2)) / 2, 1, "AbsTol", 1e-12);
            % Both triangles face the same way (+z) as the square.
            testCase.verifyTrue(all(Normals(:, 3) > 0));
        end

        function QuadsToTrisPassesTrianglesThrough(testCase)
            Verts = [0 0 0; 1 0 0; 0 1 0];
            [TriFaces, ~] = QuadsToTris([1 2 3], Verts);

            testCase.verifyEqual(TriFaces, [1 2 3]);
        end

        function QuadsToTrisRejectsBadFaceWidth(testCase)
            Verts = rand(5, 3);
            testCase.verifyError(@() QuadsToTris([1 2 3 4 5], Verts), ...
                "QuadsToTris:InvalidFaceWidth");
        end

        function QuadsToTrisRejectsOutOfRangeIndex(testCase)
            Verts = rand(3, 3);
            testCase.verifyError(@() QuadsToTris([1 2 3 4], Verts), ...
                "QuadsToTris:InvalidFaceIndex");
        end
    end

    %% Volume and centroid of known shapes

    methods (Test)
        function BoxVolumeAndCentroid(testCase)
            Corner = [1 2 3];
            BoxSize = [2 3 4];
            [Faces, Verts] = MeshUtilTests.MakeBox(Corner, BoxSize);

            testCase.verifyEqual(ComputeMeshVolume(Faces, Verts), prod(BoxSize), ...
                "RelTol", 1e-12);
            [~, Centroid] = ComputeMeshCentroid(Faces, Verts);
            testCase.verifyEqual(Centroid, Corner + BoxSize / 2, "AbsTol", 1e-12);
        end

        function CylinderVolumeAndCentroid(testCase)
            BaseCenter = [1 -2 3];
            Radius = 1.5;
            Height = 4;
            SegmentCount = 96;
            [Faces, Verts] = MeshUtilTests.MakeCylinder(BaseCenter, Radius, Height, SegmentCount);

            % Exact answer for the 96-sided prism the mesh actually is.
            PolygonArea = SegmentCount / 2 * Radius ^ 2 * sin(2 * pi / SegmentCount);
            ExpectedCentroid = BaseCenter + [0 0 Height / 2];

            Volume = ComputeMeshVolume(Faces, Verts);
            testCase.verifyEqual(Volume, PolygonArea * Height, "RelTol", 1e-12);
            [~, Centroid] = ComputeMeshCentroid(Faces, Verts);
            testCase.verifyEqual(Centroid, ExpectedCentroid, "AbsTol", 1e-10);

            % Loose check against the true cylinder. The polygon is
            % slightly smaller than the circle (about 0.07% here).
            testCase.verifyEqual(Volume, pi * Radius ^ 2 * Height, "RelTol", 2e-3);
        end

        function SphereVolumeAndCentroid(testCase)
            Center = [1 -2 3];
            Radius = 2;
            [Faces, Verts] = MeshUtilTests.MakeSphere(Center, Radius, 2000);

            % The mesh is a polyhedron inscribed in the sphere, so its
            % volume is a little under the true value.
            testCase.verifyEqual(ComputeMeshVolume(Faces, Verts), ...
                4 / 3 * pi * Radius ^ 3, "RelTol", 1e-2);
            [~, Centroid] = ComputeMeshCentroid(Faces, Verts);
            testCase.verifyEqual(Centroid, Center, "AbsTol", 0.01 * Radius);
        end

        function LPrismVolumeAndCentroid(testCase)
            Height = 2;
            [Faces, Verts] = MeshUtilTests.MakeLPrism(Height);

            testCase.verifyEqual(ComputeMeshVolume(Faces, Verts), 3 * Height, ...
                "RelTol", 1e-12);
            % The L outline is three unit squares: centroid (5/6, 5/6).
            [~, Centroid] = ComputeMeshCentroid(Faces, Verts);
            testCase.verifyEqual(Centroid, [5/6, 5/6, Height / 2], "AbsTol", 1e-12);
        end

        function CentroidIsRelativeToReferencePoint(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([1 2 3], [2 2 2]);   % center (2, 3, 4)
            ReferencePoint = [10 20 30];

            [Relative, Absolute] = ComputeMeshCentroid(Faces, Verts, ReferencePoint);

            testCase.verifyEqual(Absolute, [2 3 4], "AbsTol", 1e-12);
            testCase.verifyEqual(Relative, [2 3 4] - ReferencePoint, "AbsTol", 1e-12);
        end

        function InvertedNormalsGiveNegativeVolume(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 2 3]);
            InvertedFaces = Faces(:, [1 3 2]);

            testCase.verifyEqual(ComputeMeshVolume(InvertedFaces, Verts), -6, "RelTol", 1e-12);
        end

        function VolumeIsIndependentOfPosition(testCase)
            % Far from the origin the answer must not lose accuracy.
            [Faces, Verts] = MeshUtilTests.MakeBox([1e6 1e6 1e6], [2 3 4]);

            testCase.verifyEqual(ComputeMeshVolume(Faces, Verts), 24, "RelTol", 1e-6);
        end
    end

    %% IsMeshWatertight

    methods (Test)
        function ClosedShapesAreWatertight(testCase)
            [BoxFaces, BoxVerts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);
            [CylFaces, CylVerts] = MeshUtilTests.MakeCylinder([0 0 0], 1, 2, 32);
            [SphFaces, SphVerts] = MeshUtilTests.MakeSphere([0 0 0], 1, 500);

            testCase.verifyTrue(IsMeshWatertight(BoxFaces, BoxVerts));
            testCase.verifyTrue(IsMeshWatertight(CylFaces, CylVerts));
            testCase.verifyTrue(IsMeshWatertight(SphFaces, SphVerts));
        end

        function MissingTriangleLeavesThreeOpenEdges(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);
            Faces(1, :) = [];

            [IsWatertight, OpenEdges, NonManifoldEdges] = IsMeshWatertight(Faces, Verts);

            testCase.verifyFalse(IsWatertight);
            testCase.verifySize(OpenEdges, [3 2]);
            testCase.verifyEmpty(NonManifoldEdges);
        end

        function ExtraTriangleCreatesNonManifoldEdge(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);
            % Vertices 1 and 2 are joined by a box edge that is already
            % used twice; a third triangle on it makes it non-manifold.
            Faces = [Faces; 1 2 8];

            [IsWatertight, ~, NonManifoldEdges] = IsMeshWatertight(Faces, Verts);

            testCase.verifyFalse(IsWatertight);
            testCase.verifyTrue(ismember([1 2], NonManifoldEdges, "rows"));
        end

        function EmptyMeshIsNotWatertight(testCase)
            testCase.verifyFalse(IsMeshWatertight(zeros(0, 3), zeros(0, 3)));
        end

        function ShowFigureDrawsOnlyWhenNotWatertight(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            FiguresBefore = numel(findall(groot, "Type", "figure"));
            IsMeshWatertight(Faces, Verts, ShowFigure=true);
            testCase.verifyEqual(numel(findall(groot, "Type", "figure")), FiguresBefore, ...
                "A watertight mesh should not open a figure.");

            IsMeshWatertight(Faces(2:end, :), Verts, ShowFigure=true);
            testCase.verifyEqual(numel(findall(groot, "Type", "figure")), FiguresBefore + 1, ...
                "An open mesh should open one figure.");
        end
    end

    %% FindDegenerateTriangles and FixDegenerateTriangles

    methods (Test)
        function CleanMeshHasNoDegenerateTriangles(testCase)
            [Faces, Verts] = MeshUtilTests.MakeSphere([0 0 0], 1, 500);

            testCase.verifyFalse(any(FindDegenerateTriangles(Faces, Verts)));
        end

        function RepeatedIndexIsFlagged(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);
            FaceCount = size(Faces, 1);
            Faces = [Faces; 1 1 2];

            [IsDegenerate, HasRepeatedIndex] = FindDegenerateTriangles(Faces, Verts);

            testCase.verifyEqual(find(IsDegenerate), FaceCount + 1);
            testCase.verifyEqual(find(HasRepeatedIndex), FaceCount + 1);
        end

        function CollinearTriangleIsFlaggedWithoutRepeatedIndex(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);
            FaceCount = size(Faces, 1);
            Verts = [Verts; 0 0 0; 1 1 1; 2 2 2];   % three points on a line
            Faces = [Faces; 9 10 11];

            [IsDegenerate, HasRepeatedIndex] = FindDegenerateTriangles(Faces, Verts);

            testCase.verifyEqual(find(IsDegenerate), FaceCount + 1);
            testCase.verifyFalse(any(HasRepeatedIndex));
        end

        function CoincidentVerticesAreFlagged(testCase)
            Verts = [0 0 0; 0 0 0; 1 0 0];   % first two vertices at the same spot
            testCase.verifyTrue(FindDegenerateTriangles([1 2 3], Verts));
        end

        function ShowFigureDrawsOnlyWhenDegenerateTrianglesExist(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            FiguresBefore = numel(findall(groot, "Type", "figure"));
            FindDegenerateTriangles(Faces, Verts, ShowFigure=true);
            testCase.verifyEqual(numel(findall(groot, "Type", "figure")), FiguresBefore);

            FindDegenerateTriangles([Faces; 1 1 2], Verts, ShowFigure=true);
            testCase.verifyEqual(numel(findall(groot, "Type", "figure")), FiguresBefore + 1);
        end

        function FixRemovesRepeatedIndexTriangle(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 2 3]);
            FaceCount = size(Faces, 1);

            [FixedFaces, FixedVerts, RemovedFaceCount] = ...
                FixDegenerateTriangles([Faces; 1 1 2], Verts);

            testCase.verifyEqual(RemovedFaceCount, 1);
            testCase.verifySize(FixedFaces, [FaceCount 3]);
            testCase.verifyTrue(IsMeshWatertight(FixedFaces, FixedVerts));
            testCase.verifyEqual(ComputeMeshVolume(FixedFaces, FixedVerts), 6, "RelTol", 1e-12);
        end

        function FixCollapsesSliverTriangle(testCase)
            % Build a watertight box with one zero-area triangle in it.
            % Vertex M sits almost on top of corner A along edge A-B. The
            % triangle [A B C] is replaced by [A M C] and [M B C], and the
            % zero-area triangle [A B M] closes the gap along edge A-B.
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 2 3]);
            CornerA = Faces(1, 1);
            CornerB = Faces(1, 2);
            CornerC = Faces(1, 3);
            NewVertex = size(Verts, 1) + 1;
            Verts(NewVertex, :) = Verts(CornerA, :) + 1e-6 * (Verts(CornerB, :) - Verts(CornerA, :));
            Faces = [Faces(2:end, :); ...
                     CornerA, NewVertex, CornerC; ...
                     NewVertex, CornerB, CornerC; ...
                     CornerA, CornerB, NewVertex];

            testCase.verifyTrue(IsMeshWatertight(Faces, Verts), "Test mesh setup should be watertight.");
            testCase.verifyTrue(any(FindDegenerateTriangles(Faces, Verts)), "Test mesh should contain a sliver.");

            [FixedFaces, FixedVerts] = FixDegenerateTriangles(Faces, Verts);

            testCase.verifyFalse(any(FindDegenerateTriangles(FixedFaces, FixedVerts)));
            testCase.verifyTrue(IsMeshWatertight(FixedFaces, FixedVerts));
            testCase.verifySize(FixedFaces, [12 3]);
            testCase.verifySize(FixedVerts, [8 3]);
            % Corner A moved by about 5e-7 of an edge, so the volume is
            % still the box volume to within a tiny error.
            testCase.verifyEqual(ComputeMeshVolume(FixedFaces, FixedVerts), 6, "RelTol", 1e-5);
        end
    end

    %% FixFaceNormals

    methods (Test)
        function FixFaceNormalsRepairsPartlyFlippedBox(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 2 3]);
            WasFlipped = false(size(Faces, 1), 1);
            WasFlipped(3:3:end) = true;      % face 1 stays correct, so it seeds the fix
            BadFaces = Faces;
            BadFaces(WasFlipped, :) = BadFaces(WasFlipped, [1 3 2]);

            [FixedFaces, FixedVerts, FlipMask] = FixFaceNormals(BadFaces, Verts);

            testCase.verifyEqual(FlipMask, WasFlipped);
            testCase.verifyEqual(ComputeMeshVolume(FixedFaces, FixedVerts), 6, "RelTol", 1e-12);
            testCase.verifyEqual(FixedVerts, Verts);
        end

        function FixFaceNormalsLeavesCorrectMeshAlone(testCase)
            [Faces, Verts] = MeshUtilTests.MakeCylinder([0 0 0], 1, 2, 24);

            [FixedFaces, ~, FlipMask] = FixFaceNormals(Faces, Verts);

            testCase.verifyFalse(any(FlipMask));
            testCase.verifyEqual(FixedFaces, Faces);
        end

        function FixFaceNormalsFlipsFullyInvertedSphere(testCase)
            [Faces, Verts] = MeshUtilTests.MakeSphere([0 0 0], 1, 500);
            InvertedFaces = Faces(:, [1 3 2]);

            [FixedFaces, FixedVerts, FlipMask] = FixFaceNormals(InvertedFaces, Verts);

            testCase.verifyTrue(all(FlipMask));
            testCase.verifyGreaterThan(ComputeMeshVolume(FixedFaces, FixedVerts), 0);
        end

        function FixFaceNormalsHandlesSeparatePieces(testCase)
            % Two boxes far apart; only the second one is inside-out.
            [FacesOne, VertsOne] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);
            [FacesTwo, VertsTwo] = MeshUtilTests.MakeBox([10 0 0], [2 2 2]);
            FacesTwo = FacesTwo(:, [1 3 2]);
            Faces = [FacesOne; FacesTwo + size(VertsOne, 1)];
            Verts = [VertsOne; VertsTwo];

            [FixedFaces, FixedVerts, FlipMask] = FixFaceNormals(Faces, Verts);

            ExpectedMask = [false(size(FacesOne, 1), 1); true(size(FacesTwo, 1), 1)];
            testCase.verifyEqual(FlipMask, ExpectedMask);
            testCase.verifyEqual(ComputeMeshVolume(FixedFaces, FixedVerts), 1 + 8, "RelTol", 1e-12);
        end

        function FixFaceNormalsRejectsOpenMesh(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            testCase.verifyError(@() FixFaceNormals(Faces(2:end, :), Verts), ...
                "FixFaceNormals:NotWatertight");
        end
    end

    %% SplitMeshByPlane

    methods (Test)
        function SplitBoxAlongAxis(testCase)
            Corner = [1 2 3];
            BoxSize = [2 3 4];
            [Faces, Verts] = MeshUtilTests.MakeBox(Corner, BoxSize);

            % Plane x = 1.5 cuts the box (x from 1 to 3) into a slab 0.5
            % thick (negative side) and a slab 1.5 thick (positive side).
            [FacesPos, VertsPos, FacesNeg, VertsNeg] = ...
                SplitMeshByPlane(Faces, Verts, [1 0 0], [1.5 0 0]);

            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, 0.5 * 3 * 4, [1.25 3.5 5], 1e-9, 1e-9);
            testCase.VerifyClosedPiece(FacesPos, VertsPos, 1.5 * 3 * 4, [2.25 3.5 5], 1e-9, 1e-9);
        end

        function SplitBoxFlippedNormalSwapsSides(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([1 2 3], [2 3 4]);

            [FacesPos, VertsPos, ~, ~] = SplitMeshByPlane(Faces, Verts, [-1 0 0], [1.5 0 0]);

            % With the normal pointing toward -x, the positive piece is the
            % thin slab on the low-x side.
            testCase.VerifyClosedPiece(FacesPos, VertsPos, 0.5 * 3 * 4, [1.25 3.5 5], 1e-9, 1e-9);
        end

        function SplitBoxObliqueThroughCenter(testCase)
            % Half extents (1, 1.5, 2.3) are chosen so that no corner lies
            % on the plane x + y + z = center sum.
            Corner = [0.5 -1 2];
            BoxSize = [2 3 4.6];
            Center = Corner + BoxSize / 2;
            Normal = [1 1 1];
            [Faces, Verts] = MeshUtilTests.MakeBox(Corner, BoxSize);

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = ...
                SplitMeshByPlane(Faces, Verts, Normal, Center);

            % A plane through the center of a box splits it into two
            % mirror-image halves (point reflection through the center).
            HalfVolume = prod(BoxSize) / 2;
            testCase.verifyTrue(IsMeshWatertight(FacesPos, VertsPos));
            testCase.verifyTrue(IsMeshWatertight(FacesNeg, VertsNeg));
            testCase.verifyEqual(ComputeMeshVolume(FacesPos, VertsPos), HalfVolume, "RelTol", 1e-9);
            testCase.verifyEqual(ComputeMeshVolume(FacesNeg, VertsNeg), HalfVolume, "RelTol", 1e-9);

            [~, CentroidPos] = ComputeMeshCentroid(FacesPos, VertsPos);
            [~, CentroidNeg] = ComputeMeshCentroid(FacesNeg, VertsNeg);
            testCase.verifyEqual((CentroidPos + CentroidNeg) / 2, Center, "AbsTol", 1e-9);
            testCase.verifyGreaterThan(dot(CentroidPos - Center, Normal), 0);
            testCase.verifyLessThan(dot(CentroidNeg - Center, Normal), 0);
        end

        function SplitBoxWithPlaneThroughVertices(testCase)
            % The plane x = y passes exactly through two vertical box
            % edges, which exercises the on-plane vertex handling.
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [2 2 3]);

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = ...
                SplitMeshByPlane(Faces, Verts, [1 -1 0], [1 1 0]);

            testCase.verifyTrue(IsMeshWatertight(FacesPos, VertsPos));
            testCase.verifyTrue(IsMeshWatertight(FacesNeg, VertsNeg));
            testCase.verifyEqual(ComputeMeshVolume(FacesPos, VertsPos), 6, "RelTol", 1e-6);
            testCase.verifyEqual(ComputeMeshVolume(FacesNeg, VertsNeg), 6, "RelTol", 1e-6);
        end

        function SplitCylinderAcrossAxis(testCase)
            BaseCenter = [1 -2 3];
            Radius = 1.5;
            Height = 4;
            SegmentCount = 96;
            [Faces, Verts] = MeshUtilTests.MakeCylinder(BaseCenter, Radius, Height, SegmentCount);
            PolygonArea = SegmentCount / 2 * Radius ^ 2 * sin(2 * pi / SegmentCount);
            CutHeight = 0.3 * Height;

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = SplitMeshByPlane( ...
                Faces, Verts, [0 0 1], BaseCenter + [0 0 CutHeight]);

            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, PolygonArea * CutHeight, ...
                BaseCenter + [0 0 CutHeight / 2], 1e-9, 1e-9);
            testCase.VerifyClosedPiece(FacesPos, VertsPos, PolygonArea * (Height - CutHeight), ...
                BaseCenter + [0 0 (CutHeight + Height) / 2], 1e-9, 1e-9);
        end

        function SplitCylinderAlongAxis(testCase)
            % A plane parallel to the axis cuts the polygon outline of the
            % cylinder in two. The expected area and centroid of each part
            % come from polyshape, which is independent of the library.
            BaseCenter = [1 -2 3];
            Radius = 1.5;
            Height = 4;
            SegmentCount = 96;
            [Faces, Verts] = MeshUtilTests.MakeCylinder(BaseCenter, Radius, Height, SegmentCount);

            CutOffset = 0.3 * Radius;     % plane x = BaseCenter(1) + CutOffset
            Angles = (0:SegmentCount - 1)' * 2 * pi / SegmentCount;
            Outline = polyshape(Radius * cos(Angles), Radius * sin(Angles));
            PositiveHalfPlane = polyshape([CutOffset 4*Radius 4*Radius CutOffset], ...
                                          [-4*Radius -4*Radius 4*Radius 4*Radius]);
            NegativeHalfPlane = polyshape([-4*Radius CutOffset CutOffset -4*Radius], ...
                                          [-4*Radius -4*Radius 4*Radius 4*Radius]);
            PositivePart = intersect(Outline, PositiveHalfPlane);
            NegativePart = intersect(Outline, NegativeHalfPlane);
            [PositiveX, PositiveY] = centroid(PositivePart);
            [NegativeX, NegativeY] = centroid(NegativePart);

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = SplitMeshByPlane( ...
                Faces, Verts, [1 0 0], BaseCenter + [CutOffset 0 0]);

            testCase.VerifyClosedPiece(FacesPos, VertsPos, area(PositivePart) * Height, ...
                BaseCenter + [PositiveX, PositiveY, Height / 2], 1e-9, 1e-9);
            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, area(NegativePart) * Height, ...
                BaseCenter + [NegativeX, NegativeY, Height / 2], 1e-9, 1e-9);
        end

        function SplitSphereMakesSphericalCap(testCase)
            Center = [1 -2 3];
            Radius = 2;
            CapHeight = 0.5 * Radius;
            [Faces, Verts] = MeshUtilTests.MakeSphere(Center, Radius, 2000);

            % Cut a cap off the top of the sphere with a horizontal plane.
            [FacesPos, VertsPos, FacesNeg, VertsNeg] = SplitMeshByPlane( ...
                Faces, Verts, [0 0 1], Center + [0 0 Radius - CapHeight]);

            % Analytic cap: volume = pi*h^2*(3R - h)/3, and its centroid is
            % 3*(2R - h)^2 / (4*(3R - h)) above the sphere's center.
            CapVolume = pi * CapHeight ^ 2 * (3 * Radius - CapHeight) / 3;
            CapOffset = 3 * (2 * Radius - CapHeight) ^ 2 / (4 * (3 * Radius - CapHeight));
            SphereVolume = 4 / 3 * pi * Radius ^ 3;
            RestVolume = SphereVolume - CapVolume;
            % The sphere's centroid is its center, so the rest balances the cap.
            RestOffset = -CapVolume * CapOffset / RestVolume;

            % Loose tolerances: the mesh is a polyhedron, not a true sphere.
            testCase.VerifyClosedPiece(FacesPos, VertsPos, CapVolume, ...
                Center + [0 0 CapOffset], 3e-2, 0.03 * Radius);
            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, RestVolume, ...
                Center + [0 0 RestOffset], 3e-2, 0.03 * Radius);

            % The two pieces must add back up to the original mesh exactly.
            VolumeSum = ComputeMeshVolume(FacesPos, VertsPos) + ComputeMeshVolume(FacesNeg, VertsNeg);
            testCase.verifyEqual(VolumeSum, ComputeMeshVolume(Faces, Verts), "RelTol", 1e-9);
        end

        function SplitSphereWithTiltedPlaneMakesHemispheres(testCase)
            Center = [1 -2 3];
            Radius = 2;
            Normal = [1 2 -2] / 3;       % unit length (1 + 4 + 4 = 9)
            [Faces, Verts] = MeshUtilTests.MakeSphere(Center, Radius, 2000);

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = ...
                SplitMeshByPlane(Faces, Verts, Normal, Center);

            % A hemisphere has volume (2/3)*pi*R^3 and its centroid lies
            % 3R/8 from the center along the axis.
            HemisphereVolume = 2 / 3 * pi * Radius ^ 3;
            testCase.VerifyClosedPiece(FacesPos, VertsPos, HemisphereVolume, ...
                Center + 3 * Radius / 8 * Normal, 3e-2, 0.03 * Radius);
            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, HemisphereVolume, ...
                Center - 3 * Radius / 8 * Normal, 3e-2, 0.03 * Radius);
        end

        function SplitNonConvexPrismThroughCap(testCase)
            % Horizontal cut through the L-shaped prism: the cap is the
            % non-convex L outline.
            Height = 2;
            CutHeight = 0.5;
            [Faces, Verts] = MeshUtilTests.MakeLPrism(Height);

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = ...
                SplitMeshByPlane(Faces, Verts, [0 0 1], [0 0 CutHeight]);

            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, 3 * CutHeight, ...
                [5/6, 5/6, CutHeight / 2], 1e-9, 1e-9);
            testCase.VerifyClosedPiece(FacesPos, VertsPos, 3 * (Height - CutHeight), ...
                [5/6, 5/6, (CutHeight + Height) / 2], 1e-9, 1e-9);
        end

        function SplitNonConvexPrismIntoSeparatePieces(testCase)
            % The plane x + y = 2.5 clips off the tip of each arm of the L
            % (two separate loops, two separate pieces on the positive
            % side) and leaves the rest of the L on the negative side.
            Height = 2;
            [Faces, Verts] = MeshUtilTests.MakeLPrism(Height);

            [FacesPos, VertsPos, FacesNeg, VertsNeg] = ...
                SplitMeshByPlane(Faces, Verts, [1 1 0], [1.25 1.25 0]);

            % Each tip is a right triangle with legs 0.5 (area 0.125).
            % Their centroids are (11/6, 5/6) and (5/6, 11/6); equal areas
            % put the combined centroid at (4/3, 4/3).
            PositiveArea = 0.25;
            NegativeArea = 3 - PositiveArea;
            NegativeCentroidXY = (3 * 5/6 - PositiveArea * 4/3) / NegativeArea;   % 26/33

            testCase.VerifyClosedPiece(FacesPos, VertsPos, PositiveArea * Height, ...
                [4/3, 4/3, Height / 2], 1e-9, 1e-9);
            testCase.VerifyClosedPiece(FacesNeg, VertsNeg, NegativeArea * Height, ...
                [NegativeCentroidXY, NegativeCentroidXY, Height / 2], 1e-9, 1e-9);
        end

        function SplitRejectsPlaneThatMissesTheMesh(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            testCase.verifyError(@() SplitMeshByPlane(Faces, Verts, [1 0 0], [100 0 0]), ...
                "SplitMeshByPlane:NoIntersection");
        end

        function SplitRejectsZeroNormal(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            testCase.verifyError(@() SplitMeshByPlane(Faces, Verts, [0 0 0], [0 0 0]), ...
                "SplitMeshByPlane:ZeroNormal");
        end

        function SplitRejectsOpenMesh(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            testCase.verifyError(@() SplitMeshByPlane(Faces(2:end, :), Verts, [1 0 0], [0.5 0 0]), ...
                "SplitMeshByPlane:NotWatertight");
        end

        function SplitRejectsInvertedNormals(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            testCase.verifyError(@() SplitMeshByPlane(Faces(:, [1 3 2]), Verts, [1 0 0], [0.5 0 0]), ...
                "SplitMeshByPlane:InvertedNormals");
        end
    end

    %% VisualizeFaceNormals

    methods (Test)
        function VisualizeFaceNormalsColorsFacesByVisibility(testCase)
            [Faces, Verts] = MeshUtilTests.MakeBox([0 0 0], [1 1 1]);

            FigureHandle = VisualizeFaceNormals(Faces, Verts);
            drawnow;
            Ax = findobj(FigureHandle, "Type", "axes");
            PatchHandle = findobj(FigureHandle, "Type", "patch");
            testCase.assumeTrue(strcmp(Ax.Projection, "orthographic"));

            % Work out independently which triangles face the camera.
            CornerA = Verts(Faces(:, 1), :);
            CornerB = Verts(Faces(:, 2), :);
            CornerC = Verts(Faces(:, 3), :);
            Normals = cross(CornerB - CornerA, CornerC - CornerA, 2);
            ExpectedFacing = Normals * (Ax.CameraPosition - Ax.CameraTarget)' > 0;

            Colors = PatchHandle.FaceVertexCData;
            IsGreen = Colors(:, 2) > Colors(:, 1);

            testCase.verifyEqual(IsGreen, ExpectedFacing);
            % From the default oblique 3-D view a box shows three of its
            % six sides, i.e. six of its twelve triangles.
            testCase.verifyEqual(nnz(IsGreen), 6);
        end
    end

    %% Helpers

    methods (Access = private)
        function VerifyClosedPiece(testCase, Faces, Verts, ExpectedVolume, ExpectedCentroid, VolumeRelTol, CentroidAbsTol)
            %VERIFYCLOSEDPIECE Check a mesh is closed, outward-facing, and has the expected volume and centroid.

            % Assert (stop this test) rather than verify, since the checks
            % below assume a closed mesh and FixFaceNormals errors otherwise.
            testCase.assertTrue(IsMeshWatertight(Faces, Verts), "Piece should be watertight.");

            % Correct normals means FixFaceNormals has nothing to flip.
            [~, ~, FlipMask] = FixFaceNormals(Faces, Verts);
            testCase.verifyFalse(any(FlipMask), "Piece should have consistent outward normals.");

            testCase.verifyEqual(ComputeMeshVolume(Faces, Verts), ExpectedVolume, ...
                "RelTol", VolumeRelTol);
            [~, Centroid] = ComputeMeshCentroid(Faces, Verts);
            testCase.verifyEqual(Centroid, ExpectedCentroid, "AbsTol", CentroidAbsTol);
        end
    end

    methods (Static, Access = private)
        function [Faces, Verts] = MakeBox(Corner, BoxSize)
            %MAKEBOX Axis-aligned box with its minimum corner at Corner. 12 triangles.

            % Corner k has x, y, z offsets given by the bits of k - 1:
            % index = 1 + x + 2*y + 4*z.
            Bits = [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];
            Verts = Corner + Bits .* BoxSize;

            % The six sides, each listed in order around the side.
            Quads = [1 3 7 5;     % x = min
                     2 4 8 6;     % x = max
                     1 2 6 5;     % y = min
                     3 4 8 7;     % y = max
                     1 2 4 3;     % z = min
                     5 6 8 7];    % z = max
            Faces = [Quads(:, [1 2 3]); Quads(:, [1 3 4])];

            Faces = MeshUtilTests.OrientOutward(Faces, Verts, Corner + BoxSize / 2);
        end

        function [Faces, Verts] = MakeCylinder(BaseCenter, Radius, Height, SegmentCount)
            %MAKECYLINDER Closed cylinder along +z with flat caps (a SegmentCount-sided prism).

            Angles = (0:SegmentCount - 1)' * 2 * pi / SegmentCount;
            RingX = BaseCenter(1) + Radius * cos(Angles);
            RingY = BaseCenter(2) + Radius * sin(Angles);
            BottomRing = [RingX, RingY, repmat(BaseCenter(3), SegmentCount, 1)];
            TopRing = [RingX, RingY, repmat(BaseCenter(3) + Height, SegmentCount, 1)];

            % Vertices: bottom ring, top ring, bottom center, top center.
            Verts = [BottomRing; TopRing; BaseCenter; BaseCenter + [0 0 Height]];
            BottomCenter = 2 * SegmentCount + 1;
            TopCenter = 2 * SegmentCount + 2;

            Current = (1:SegmentCount)';
            Next = mod(Current, SegmentCount) + 1;

            SideFaces = [Current, Next, SegmentCount + Next; ...
                         Current, SegmentCount + Next, SegmentCount + Current];
            TopFaces = [repmat(TopCenter, SegmentCount, 1), SegmentCount + Current, SegmentCount + Next];
            BottomFaces = [repmat(BottomCenter, SegmentCount, 1), Next, Current];
            Faces = [SideFaces; TopFaces; BottomFaces];

            Faces = MeshUtilTests.OrientOutward(Faces, Verts, BaseCenter + [0 0 Height / 2]);
        end

        function [Faces, Verts] = MakeSphere(Center, Radius, PointCount)
            %MAKESPHERE Polyhedron through PointCount evenly spread points on a sphere.

            % Fibonacci spiral: evenly spaced heights, golden-angle turns.
            Index = (0:PointCount - 1)';
            Z = 1 - 2 * (Index + 0.5) / PointCount;
            RingRadius = sqrt(1 - Z .^ 2);
            Phi = Index * pi * (3 - sqrt(5));
            Verts = Center + Radius * [RingRadius .* cos(Phi), RingRadius .* sin(Phi), Z];

            % Every point is on a sphere, so all are corners of the convex
            % hull and the hull's triangles form the closed surface.
            Faces = convhull(Verts(:, 1), Verts(:, 2), Verts(:, 3));

            Faces = MeshUtilTests.OrientOutward(Faces, Verts, Center);
        end

        function [Faces, Verts] = MakeLPrism(Height)
            %MAKELPRISM L-shaped prism (non-convex), outline area 3, standing on z = 0.

            % Counter-clockwise L outline made of three unit squares.
            Outline = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
            Verts = [Outline, zeros(6, 1); Outline, Height * ones(6, 1)];

            % Fan of triangles from outline corner 1 covers the L exactly.
            OutlineTris = [1 2 3; 1 3 4; 1 4 5; 1 5 6];
            TopFaces = OutlineTris + 6;               % counter-clockwise: faces +z
            BottomFaces = OutlineTris(:, [1 3 2]);    % reversed: faces -z

            % One rectangle (two triangles) per outline edge, wound so the
            % normal points away from the inside of the L.
            Current = (1:6)';
            Next = mod(Current, 6) + 1;
            SideFaces = [Current, Next, 6 + Next; ...
                         Current, 6 + Next, 6 + Current];

            Faces = [TopFaces; BottomFaces; SideFaces];
        end

        function Faces = OrientOutward(Faces, Verts, InteriorPoint)
            %ORIENTOUTWARD Flip triangles so their normals point away from a point inside a convex shape.

            CornerA = Verts(Faces(:, 1), :);
            CornerB = Verts(Faces(:, 2), :);
            CornerC = Verts(Faces(:, 3), :);
            Normals = cross(CornerB - CornerA, CornerC - CornerA, 2);
            FaceCenters = (CornerA + CornerB + CornerC) / 3;

            PointsInward = dot(Normals, FaceCenters - InteriorPoint, 2) < 0;
            Faces(PointsInward, :) = Faces(PointsInward, [1 3 2]);
        end
    end
end
