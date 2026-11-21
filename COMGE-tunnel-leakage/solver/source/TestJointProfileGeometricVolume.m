classdef TestJointProfileGeometricVolume < matlab.unittest.TestCase
    % Exact radial geometry must not depend on midpoint thermal quadrature.
    methods (Test)
        function positiveSegmentsUseTrapezoidsAndOriginalNodePartials(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;1.5;2],[1 2 3]);
            [V,d,D]=joint_profile_geometric_volume(B,u,2e-10,.2);
            testCase.verifyEqual(V,.4,'AbsTol',1e-15);
            testCase.verifyEqual(d,[.05 .1 .05],'AbsTol',1e-15);
            testCase.verifyEqual(D.segmentVolume_m2,[.15 .25],'AbsTol',1e-15);
            testCase.verifyEqual(D.nodalAperture_m,[.2 .4 .6],'AbsTol',1e-15);
        end

        function whollyNonpositiveProfileHasNoGeometricInventory(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;1.5;2],[-1 -2 0]);
            [V,d,D,du]=joint_profile_geometric_volume(B,u,2e-10,.2);
            testCase.verifyEqual(V,0,'AbsTol',0);
            testCase.verifyEqual(d,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(D.nodalAperture_m,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(full(du),zeros(1,6),'AbsTol',0);
        end

        function positiveToNegativeCrossingHasExactTriangularArea(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;1.6],[2 -1]);
            [V,d]=joint_profile_geometric_volume(B,u,2e-10,.25);
            testCase.verifyEqual(V,.1,'AbsTol',2e-15);
            testCase.verifyEqual(d,[1/15 1/30],'AbsTol',2e-15);
        end

        function negativeToPositiveCrossingReversesNodePartials(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;1.6],[-1 2]);
            [V,d]=joint_profile_geometric_volume(B,u,2e-10,.25);
            testCase.verifyEqual(V,.1,'AbsTol',2e-15);
            testCase.verifyEqual(d,[1/30 1/15],'AbsTol',2e-15);
        end

        function zeroEndpointNextToPositiveEndpointHasContinuousDerivative(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;2],[0 2]);
            [V,d]=joint_profile_geometric_volume(B,u,2e-10,.1);
            testCase.verifyEqual(V,.1,'AbsTol',1e-15);
            testCase.verifyEqual(d,[.05 .05],'AbsTol',1e-15);
        end

        function tinyPositiveEndpointMissedByEightThermalMidpointsIsRetained(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;2],[1e-4 -.01]);
            samples=2-((1:8)'-.5)/8;
            midpointGap=1e-4+(-.01-1e-4)*(samples-1);
            [V,~,D]=joint_profile_geometric_volume(B,u,2e-10,1);
            testCase.verifyLessThan(midpointGap,zeros(8,1));
            testCase.verifyEqual(V,.5*1e-8/.0101,'RelTol',2e-14);
            testCase.verifyGreaterThan(V,0);
            testCase.verifyFalse(D.wholePathOpen);
        end

        function descendingRadialOrderPreservesExactAreaAndMappedDerivatives(testCase)
            [A,ua]=TestJointProfileGeometricVolume.fixture([1;1.4;2],[.002 -.001 .004]);
            [B,ub]=TestJointProfileGeometricVolume.fixture([2;1.4;1],[.004 -.001 .002]);
            [va,da]=joint_profile_geometric_volume(A,ua,2e-10,.001);
            [vb,db]=joint_profile_geometric_volume(B,ub,2e-10,.001);
            testCase.verifyEqual(vb,va,'AbsTol',1e-20);
            testCase.verifyEqual(db,fliplr(da),'AbsTol',1e-18);
        end

        function nodalGapPartialsMatchCentralFiniteDifferences(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([3;2.4;1],[.002 -.001 .004]);
            [~,analytic]=joint_profile_geometric_volume(B,u,2e-10,.001);
            finite=TestJointProfileGeometricVolume.finiteGapJacobian(B,u,2e-10,.001,1e-8);
            testCase.verifyEqual(analytic,finite,'AbsTol',2e-13,'RelTol',2e-9);
        end

        function ringDisplacementDerivativeUsesJointMajorGapMap(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([3;2.4;1],[.002 -.001 .004]);
            direction=[.1;-.3;.2;0;0;0]; epsilon=1e-8;
            [~,d,~,du]=joint_profile_geometric_volume(B,u,2e-10,.001);
            plus=joint_profile_geometric_volume(B,u+epsilon*direction,2e-10,.001);
            minus=joint_profile_geometric_volume(B,u-epsilon*direction,2e-10,.001);
            testCase.verifyTrue(issparse(du));
            testCase.verifyEqual(full(du),[d zeros(1,3)],'AbsTol',0);
            testCase.verifyEqual(full(du*direction),(plus-minus)/(2*epsilon), ...
                'AbsTol',2e-13,'RelTol',2e-9);
        end

        function uniformProfileUsesActualPerJointPathLengths(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture( ...
                [3 4;2.5 3.5;1 3],[.002 .002 .002;.004 .004 .004]);
            [V,d,D]=joint_profile_geometric_volume(B,u,[2e-10;8e-10],[.1;.2]);
            testCase.verifyEqual(V,[.0004;.0008],'AbsTol',1e-18);
            testCase.verifyEqual(D.actualPathLength_m,[2;1],'AbsTol',0);
            testCase.verifyEqual(d,[.025 .1 .075;.05 .1 .05],'AbsTol',1e-16);
        end

        function sealedAndDisabledJointsHaveExactlyZeroVolume(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture( ...
                [3 3;2 2;1 1],[.002 -.001 .004;.003 .003 .003]);
            [V,d,D,du]=joint_profile_geometric_volume(B,u,[0;2e-10],[1;0]);
            testCase.verifyEqual(V,zeros(2,1),'AbsTol',0);
            testCase.verifyEqual(d,zeros(2,3),'AbsTol',0);
            testCase.verifyEqual(D.nodalAperture_m,zeros(2,3),'AbsTol',0);
            testCase.verifyEqual(full(du),zeros(2,12),'AbsTol',0);
            testCase.verifyEqual(D.sealedMask,[true;false]);
        end

        function realGapProfileCacheAndUniformQuadratureVolumeAgree(testCase)
            S=make_thm_verification_mesh(16); layout=hsr_eight_segment_layout();
            M=build_segment_mechanical_mesh(S,layout.jointAnglesRad);
            C=solve_segment_contact_fem('assemble',M,struct());
            C.gapNodal0_m=.002*ones(numel(C.gapNodal0_m),1);
            paths=build_joint_thermal_paths(S,layout.jointAnglesRad,8);
            B=build_joint_gap_profiles(C,paths); u=zeros(size(B.nodalGapU,2),1);
            [V,~,D]=joint_profile_geometric_volume(B,u,2e-10,.001);
            options=struct('pathLength_m',B.actualPathLength, ...
                'openingTransmissionFraction',.001,'sampleWeights',B.sampleWeights);
            [~,~,~,thermal]=joint_opening_profile_hydraulics(2e-10,B.gapProfile_m,[],options);
            testCase.verifyEqual(V,B.actualPathLength*2e-6,'AbsTol',1e-19);
            testCase.verifyEqual(V,thermal.newslitGeometricVolume_m2,'AbsTol',1e-19);
            testCase.verifyEqual(D.actualPathLength_m,B.actualPathLength,'AbsTol',1e-14);
        end

        function physicalParameterValidationIsIndependentOfSealedMask(testCase)
            [B,u]=TestJointProfileGeometricVolume.fixture([1;2],[1 1]);
            testCase.verifyError(@()joint_profile_geometric_volume(B,u,-1,.1), ...
                'JointProfileVolume:Conductance');
            testCase.verifyError(@()joint_profile_geometric_volume(B,u,0,1.1), ...
                'JointProfileVolume:Fraction');
            B.nodalRadius_m=[1;1];
            testCase.verifyError(@()joint_profile_geometric_volume(B,u,0,.1), ...
                'JointProfileVolume:Radius');
        end
    end

    methods (Static,Access=private)
        function [B,u]=fixture(radius,gap)
            nr=size(radius,1); nj=size(gap,1);
            B=struct('nr',nr,'nj',nj,'nodalRadius_m',radius, ...
                'nodalGap0',reshape(gap',[],1), ...
                'nodalGapU',sparse((1:nr*nj)',(1:nr*nj)',1,nr*nj,2*nr*nj));
            u=zeros(2*nr*nj,1);
        end

        function J=finiteGapJacobian(B,u,c,f,epsilon)
            J=zeros(B.nj,B.nr);
            for k=1:B.nr
                plus=B; minus=B;
                plus.nodalGap0(k)=plus.nodalGap0(k)+epsilon;
                minus.nodalGap0(k)=minus.nodalGap0(k)-epsilon;
                J(:,k)=(joint_profile_geometric_volume(plus,u,c,f) ...
                    -joint_profile_geometric_volume(minus,u,c,f))/(2*epsilon);
            end
        end
    end
end
