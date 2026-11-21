classdef TestSoilSmallStrainDiagnostics < matlab.unittest.TestCase
    % Complete T3 diagnostics expose shear and rotation beyond div(u).
    methods (Test)
        function affineTensionCompressionAndPrincipalStrains(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            u=[.02*xy(:,1)+.7,-.03*xy(:,2)-.9];
            [d,A]=soil_small_strain_diagnostics(xy,tri,u);
            testCase.verifyEqual(d.strainComponents,repmat([.02 -.03 0],4,1),'AbsTol',1e-14);
            testCase.verifyEqual(d.principalStrains,repmat([-.03 .02 0],4,1),'AbsTol',1e-14);
            testCase.verifyEqual(d.maxAbsVolumetricStrain,.01,'AbsTol',1e-14);
            testCase.verifyEqual(d.maxAbsPrincipalStrain,.03,'AbsTol',1e-14);
            testCase.verifyEqual(d.maximumTensorShearStrain,.025,'AbsTol',1e-14);
            testCase.verifyEqual(d.maxAbsTensorShearXY,0,'AbsTol',1e-14);
            testCase.verifyEqual(A.V*reshape(u',[],1),-.01*ones(4,1),'AbsTol',1e-14);
            testCase.verifyTrue(issparse(A.B)&&issparse(A.gradient)&&issparse(A.V));
        end

        function pureSymmetricShearHasZeroVolumeAndNonzeroPrincipalStrain(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            u=[.04*xy(:,2),.04*xy(:,1)];
            d=soil_small_strain_diagnostics(xy,tri,u);
            testCase.verifyEqual(d.strainComponents,repmat([0 0 .08],4,1),'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsVolumetricStrain,0,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsPrincipalStrain,.04,'AbsTol',1e-15);
            testCase.verifyEqual(d.maximumTensorShearStrain,.04,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsTensorShearXY,.04,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxRotationNorm,0,'AbsTol',1e-15);
        end

        function simpleShearSeparatesTensorShearAndRotation(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            u=[.08*xy(:,2),zeros(size(xy,1),1)];
            d=soil_small_strain_diagnostics(xy,tri,u);
            testCase.verifyEqual(d.tensorShearXY,.04*ones(4,1),'AbsTol',1e-15);
            testCase.verifyEqual(d.rotationComponent,-.04*ones(4,1),'AbsTol',1e-15);
            testCase.verifyEqual(d.maxGradientNorm,.08,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxRotationNorm,sqrt(2)*.04,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsVolumetricStrain,0,'AbsTol',1e-15);
        end

        function infinitesimalRigidRotationHasNoStrain(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            omega=.12; u=[-omega*xy(:,2)+3,omega*xy(:,1)-5];
            d=soil_small_strain_diagnostics(xy,tri,u);
            testCase.verifyEqual(d.strainComponents,zeros(4,3),'AbsTol',3e-14);
            testCase.verifyEqual(d.maxAbsPrincipalStrain,0,'AbsTol',3e-14);
            testCase.verifyEqual(d.maxAbsVolumetricStrain,0,'AbsTol',3e-14);
            testCase.verifyEqual(d.rotationComponent,omega*ones(4,1),'AbsTol',3e-14);
            testCase.verifyEqual(d.maxGradientNorm,sqrt(2)*omega,'AbsTol',3e-14);
            testCase.verifyEqual(d.maxRotationNorm,sqrt(2)*omega,'AbsTol',3e-14);
        end

        function rigidTranslationHasNoGradient(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            d=soil_small_strain_diagnostics(xy,tri,repmat([.8 -1.3],size(xy,1),1));
            testCase.verifyEqual(d.displacementGradient,zeros(4,4),'AbsTol',1e-14);
            testCase.verifyEqual(d.maxGradientNorm,0,'AbsTol',1e-14);
            testCase.verifyEqual(d.maxRotationNorm,0,'AbsTol',1e-14);
        end

        function planeStrainDilationIncludesOutOfPlaneZero(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            d=soil_small_strain_diagnostics(xy,tri,.01*xy);
            testCase.verifyEqual(d.principalStrains,repmat([.01 .01 0],4,1),'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsPrincipalStrain,.01,'AbsTol',1e-15);
            testCase.verifyEqual(d.maximumTensorShearStrain,.005,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsTensorShearXY,0,'AbsTol',1e-15);
            testCase.verifyEqual(d.maxAbsVolumetricStrain,.02,'AbsTol',1e-15);
        end

        function reversedTriangleOrientationIsInvariant(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            u=[.02*xy(:,1)-.03*xy(:,2),.04*xy(:,1)+.01*xy(:,2)];
            [a,A]=soil_small_strain_diagnostics(xy,tri,u);
            [b,B]=soil_small_strain_diagnostics(xy,tri(:,[1 3 2]),u);
            testCase.verifyEqual(a.strainComponents,b.strainComponents,'AbsTol',1e-15);
            testCase.verifyEqual(a.displacementGradient,b.displacementGradient,'AbsTol',1e-15);
            testCase.verifyEqual(A.B,B.B,'AbsTol',1e-15);
            testCase.verifyEqual(A.gradient,B.gradient,'AbsTol',1e-15);
            testCase.verifyEqual(A.elementAreas,B.elementAreas,'AbsTol',1e-15);
        end

        function cachedEvaluationMatchesFreshGeometry(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            [~,A]=soil_small_strain_diagnostics(xy,tri,zeros(2*size(xy,1),1));
            u=[.02*xy(:,1)+.03*xy(:,2),-.01*xy(:,1)-.04*xy(:,2)];
            fresh=soil_small_strain_diagnostics(xy,tri,u);
            cached=soil_small_strain_diagnostics(A,u);
            testCase.verifyEqual(cached,fresh);
        end

        function interleavedAndMatrixInputsAreEquivalent(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            u=[-.02*xy(:,1)+.03*xy(:,2),.06*xy(:,1)-.04*xy(:,2)];
            matrix=soil_small_strain_diagnostics(xy,tri,u);
            vector=soil_small_strain_diagnostics(xy,tri,reshape(u',1,[]));
            testCase.verifyEqual(matrix,vector);
        end

        function degenerateTriangleIsRejected(testCase)
            testCase.verifyError(@()soil_small_strain_diagnostics( ...
                [0 0;1 0;2 0],[1 2 3],zeros(6,1)),'SoilSmallStrain:Degenerate');
        end

        function wrongDisplacementSizeIsRejected(testCase)
            [xy,tri]=TestSoilSmallStrainDiagnostics.fixture();
            testCase.verifyError(@()soil_small_strain_diagnostics(xy,tri,zeros(9,1)), ...
                'SoilSmallStrain:Displacement');
        end
    end
    methods (Static,Access=private)
        function [xy,tri]=fixture()
            xy=[0 0;1 0;0 1;1 1;.5 .5];
            tri=[1 2 5;2 4 5;4 3 5;3 1 5];
        end
    end
end
