classdef TestJointGapProfiles < matlab.unittest.TestCase
    %TESTJOINTGAPPROFILES Public geometry interpolation and closure evidence.
    % Synthetic coincident faces retain separate displacement DOFs.

    methods (Test)
        function constantGapRemainsConstantWithAscendingMechanicalRadii(testCase)
            [C,paths]=TestJointGapProfiles.fixture([1;2;3],[pi/2;0;pi],4);
            C.gapNodal0_m=.002*ones(9,1);
            B=build_joint_gap_profiles(C,paths);
            testCase.verifyEqual(B.gapProfile_m,.002*ones(3,4),'AbsTol',1e-16);
            testCase.verifyEqual(full(sum(B.P,2)),ones(12,1),'AbsTol',1e-14);
            testCase.verifyGreaterThanOrEqual(nonzeros(B.P),0);
            testCase.verifyEqual(B.actualPathLength,2*ones(3,1),'AbsTol',1e-14);
            testCase.verifyEqual(B.sampleWeights,.25*ones(1,4),'AbsTol',0);
            testCase.verifyTrue(all(B.wholePathOpen));
        end

        function linearRadialGapIsExactAndDisplacementMapAgrees(testCase)
            [C,paths,U]=TestJointGapProfiles.linearFixture();
            B=build_joint_gap_profiles(C,paths,U);
            expected=.0002+.0007*B.sampleRadius_m';
            nodal=reshape(C.gapNodal0_m+C.gapNodalMap*U,3,3);
            testCase.verifyEqual(B.gapProfile_m,expected,'AbsTol',1e-15);
            testCase.verifyEqual(B.profileGapU*U,B.P*C.gapNodalMap*U,'AbsTol',1e-15);
            testCase.verifyEqual(B.nodalGapEvidence,nodal,'AbsTol',0);
            testCase.verifyEqual(B.P*B.nodalXY,paths.sampleXY,'AbsTol',1e-14);
            testCase.verifyTrue(all(B.wholePathOpen));
        end

        function duplicateCoordinatePathsRetainIndependentUnorderedJointIDs(testCase)
            [C,paths]=TestJointGapProfiles.fixture([3;2;1],[pi;0;pi],4);
            C.gapNodal0_m=repelem([.001;.002;.003],3);
            B=build_joint_gap_profiles(C,paths);
            testCase.verifyEqual(B.gapProfile_m,repmat([.001;.002;.003],1,4),'AbsTol',1e-16);
            testCase.verifyEqual(full(B.P(1:4,4:9)),zeros(4,6),'AbsTol',0);
            testCase.verifyEqual(full(B.P(5:8,[1:3 7:9])),zeros(4,6),'AbsTol',0);
            testCase.verifyEqual(paths.sampleXY(1:4,:),paths.sampleXY(9:12,:),'AbsTol',0);
        end

        function endpointsUseExactNodalValuesWithoutExtrapolation(testCase)
            [C,paths]=TestJointGapProfiles.fixture([3;2;1],0,3);
            paths.sampleXY=[4 -2;3 -2;2 -2]; % centre [1,-2], radii 3,2,1
            C.gapNodal0_m=[.001;.003;.007];
            B=build_joint_gap_profiles(C,paths);
            testCase.verifyEqual(full(B.P),eye(3),'AbsTol',0);
            testCase.verifyEqual(B.gapProfile_m,[.001 .003 .007],'AbsTol',0);
        end

        function closedEndpointIsDetectedEvenWhenEveryMidpointIsPositive(testCase)
            [C,paths]=TestJointGapProfiles.fixture([3;2;1],0,2);
            C.gapNodal0_m=[0;.002;.002];
            B=build_joint_gap_profiles(C,paths);
            testCase.verifyGreaterThan(B.gapProfile_m,zeros(1,2));
            testCase.verifyFalse(B.wholePathOpen);
            testCase.verifyTrue(B.samplingMissesClosure);
            testCase.verifyEqual(B.nodalGapEvidence(1),0,'AbsTol',0);
        end

        function closedInteriorNodeIsDetectedWithoutDestroyingPositiveInventory(testCase)
            [C,paths]=TestJointGapProfiles.fixture([3;2;1],0,2);
            C.gapNodal0_m=[.002;0;.002];
            B=build_joint_gap_profiles(C,paths);
            options=struct('openingTransmissionFraction',.001,'openPathMask',B.wholePathOpen, ...
                'pathLength_m',B.actualPathLength,'sampleWeights',B.sampleWeights);
            [value,~,~,D]=joint_opening_profile_hydraulics(2e-10,B.gapProfile_m,[2 3],options);
            testCase.verifyFalse(B.wholePathOpen);
            testCase.verifyGreaterThan(B.gapProfile_m,zeros(1,2));
            testCase.verifyEqual(value,2e-10,'AbsTol',0);
            testCase.verifyEqual(D.openConductance,0,'AbsTol',0);
            testCase.verifyGreaterThan(D.newslitGeometricVolume_m2,0);
        end

        function outOfDomainAndInconsistentEndpointAreRejected(testCase)
            [C,paths]=TestJointGapProfiles.fixture([3;2;1],0,2);
            outside=paths; outside.sampleXY(1,1)=4.1;
            mismatch=paths; mismatch.radiusOuter_m=3.1;
            testCase.verifyError(@()build_joint_gap_profiles(C,outside),'JointGapProfiles:OutsidePath');
            testCase.verifyError(@()build_joint_gap_profiles(C,mismatch),'JointGapProfiles:PathBoundary');
        end

        function mismatchedRayAndReversedSampleOrderAreRejected(testCase)
            [C,paths]=TestJointGapProfiles.fixture([3;2;1],0,2);
            offRay=paths; offRay.sampleXY(1,2)=offRay.sampleXY(1,2)+.01;
            reverse=paths; reverse.sampleXY=flipud(reverse.sampleXY);
            testCase.verifyError(@()build_joint_gap_profiles(C,offRay),'JointGapProfiles:Ray');
            testCase.verifyError(@()build_joint_gap_profiles(C,reverse),'JointGapProfiles:SampleOrder');
        end

        function realPolygonalMechanicalAndThermalPathsCoincide(testCase)
            source=make_thm_verification_mesh(16);
            layout=hsr_eight_segment_layout(); order=[7 2 5 1 8 3 6 4];
            angles=layout.jointAnglesRad(order);
            M=build_segment_mechanical_mesh(source,angles);
            C=solve_segment_contact_fem('assemble',M,struct());
            paths=build_joint_thermal_paths(source,angles,8);
            B=build_joint_gap_profiles(C,paths);
            testCase.verifyEqual(B.P*B.nodalXY,paths.sampleXY,'AbsTol',1e-9);
            testCase.verifyEqual(B.actualPathLength,paths.radiusOuter_m-paths.radiusInner_m,'AbsTol',0);
            testCase.verifySize(B.profileGapU,[64 2*size(M.p,1)]);
            testCase.verifyEqual(full(sum(B.P,2)),ones(64,1),'AbsTol',1e-12);
        end
    end

    methods (Static,Access=private)
        function [C,paths]=fixture(radii,angles,nq)
            radii=radii(:); angles=angles(:); nr=numel(radii); nj=numel(angles);
            center=[1 -2]; xy=zeros(2*nr*nj,2); pairs=zeros(nr,2,nj);
            row=zeros(4*nr*nj,1); column=row; value=row; cursor=0;
            samples=zeros(nj*nq,2);
            for j=1:nj
                dir=[cos(angles(j)),sin(angles(j))]; normal=[-sin(angles(j)),cos(angles(j))];
                minus=2*nr*(j-1)+(1:nr)'; plus=minus+nr;
                xy(minus,:)=center+radii*dir; xy(plus,:)=xy(minus,:);
                pairs(:,:,j)=[minus plus];
                for k=1:nr
                    id=(j-1)*nr+k; dofs=[2*minus(k)-1;2*minus(k);2*plus(k)-1;2*plus(k)];
                    slots=cursor+(1:4); row(slots)=id; column(slots)=dofs;
                    value(slots)=[-normal';normal']; cursor=cursor+4;
                end
                outer=max(radii); inner=min(radii);
                query=outer-((1:nq)'-.5)/nq*(outer-inner);
                samples((j-1)*nq+(1:nq),:)=center+query*dir;
            end
            M=struct('p',xy,'center',center,'jointNodePairs',pairs,'jointAnglesRad',angles);
            C=struct('mesh',M,'gapNodalMap',sparse(row,column,value,nr*nj,2*size(xy,1)), ...
                'gapNodal0_m',zeros(nr*nj,1));
            paths=struct('nj',nj,'nq',nq,'sampleXY',samples, ...
                'radiusOuter_m',max(radii)*ones(nj,1),'radiusInner_m',min(radii)*ones(nj,1));
        end

        function [C,paths,U]=linearFixture()
            radii=[3;2;1]; angles=[pi/2;0;pi];
            [C,paths]=TestJointGapProfiles.fixture(radii,angles,4);
            C.gapNodal0_m=.0002*ones(9,1); U=zeros(size(C.gapNodalMap,2),1);
            for j=1:3
                normal=[-sin(angles(j));cos(angles(j))];
                plus=C.mesh.jointNodePairs(:,2,j);
                for k=1:3
                    dofs=[2*plus(k)-1;2*plus(k)]; U(dofs)=.0007*radii(k)*normal;
                end
            end
        end
    end
end
