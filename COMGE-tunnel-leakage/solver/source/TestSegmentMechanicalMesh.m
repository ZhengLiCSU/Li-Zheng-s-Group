classdef TestSegmentMechanicalMesh < matlab.unittest.TestCase
    %TESTSEGMENTMECHANICALMESH Independent eight-block geometry checks.
    % Run results=runtests('TestSegmentMechanicalMesh.m'); assertSuccess(results).
    % No contact law or mechanical solution is asserted by these mesh tests.

    properties (Access=private)
        Source
        Angles
        Mechanical
        Report
    end

    methods (TestClassSetup)
        function buildActualSourceMesh(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.Source=load(fullfile(root,'example_mesh_repaired_interface.mat'));
            layout=hsr_eight_segment_layout();
            testCase.Angles=layout.jointAnglesRad;
            [testCase.Mechanical,testCase.Report]=build_segment_mechanical_mesh( ...
                testCase.Source,testCase.Angles);
        end
    end

    methods (Test)
        function eightDisconnectedMechanicalBlocks(testCase)
            M=testCase.Mechanical;
            components=TestSegmentMechanicalMesh.graphComponents(M.q,size(M.p,1));
            testCase.verifyEqual(numel(unique(components)),8);
            testCase.verifyEqual(numel(unique(M.nodeSegmentID)),8);
            testCase.verifyEqual(testCase.Report.disconnectedComponents,8);
        end

        function noQ4BridgesTwoSegments(testCase)
            M=testCase.Mechanical;
            nodeLabels=M.nodeSegmentID(M.q);
            testCase.verifyEqual(nodeLabels,repmat(M.elementSegmentID,1,4));
            testCase.verifyEqual(numel(unique(M.q(:))),size(M.p,1));
        end

        function coincidentJointFacesHaveIndependentDOFs(testCase)
            M=testCase.Mechanical;
            minus=reshape(M.jointNodePairs(:,1,:),[],1);
            plus=reshape(M.jointNodePairs(:,2,:),[],1);
            testCase.verifyTrue(all(minus~=plus));
            testCase.verifyEqual(M.p(minus,:),M.p(plus,:),'AbsTol',1e-10);
            testCase.verifyTrue(all(M.nodeSegmentID(minus)~=M.nodeSegmentID(plus)));
            testCase.verifyEqual(size(M.jointNodePairs,3),8);
        end

        function affineThermalPullbackIsExact(testCase)
            M=testCase.Mechanical; S=testCase.Source;
            oldTemperature=7+.31*S.pAll(:,1)-.43*S.pAll(:,2);
            expected=7+.31*M.p(:,1)-.43*M.p(:,2);
            testCase.verifyEqual(M.thermalInterpolation*oldTemperature,expected,'AbsTol',1e-11);
            testCase.verifyEqual(full(sum(M.thermalInterpolation,2)), ...
                ones(size(M.p,1),1),'AbsTol',1e-13);
            testCase.verifyGreaterThanOrEqual(nonzeros(M.thermalInterpolation),0);
            testCase.verifySize(M.thermalInterpolation,[size(M.p,1),size(S.pAll,1)]);
        end

        function originalQ4FieldRestrictionIsExact(testCase)
            M=testCase.Mechanical; S=testCase.Source;
            temperature=sin((1:size(S.pAll,1))'*.137);
            [child,source]=TestSegmentMechanicalMesh.quadraturePullback(M,S,temperature);
            testCase.verifyEqual(child,source,'AbsTol',2e-10);
        end

        function jointDirectionsAndRadialOrderAreConsistent(testCase)
            M=testCase.Mechanical; I=M.jointInterface;
            angles=M.jointAnglesRad(I.jointID);
            expectedNormal=[-sin(angles),cos(angles)];
            expectedTangent=[-cos(angles),-sin(angles)];
            testCase.verifyEqual(I.normalMinusToPlus,expectedNormal,'AbsTol',1e-13);
            testCase.verifyEqual(I.tangentOuterToInner,expectedTangent,'AbsTol',1e-10);
            testCase.verifyGreaterThan(I.length_m,0);
            testCase.verifyEqual(I.segmentPairs,M.jointSegmentPairs(I.jointID,:));
        end

        function positiveElementJacobianAndConservedArea(testCase)
            M=testCase.Mechanical;
            minimum=TestSegmentMechanicalMesh.minimumQuadJacobian(M.p,M.q);
            testCase.verifyGreaterThan(minimum,0);
            testCase.verifyEqual(testCase.Report.mechanicalShellArea_m2, ...
                testCase.Report.sourceShellArea_m2,'AbsTol',1e-9);
        end

        function soilInterfacePullbackRetainsExistingSoil(testCase)
            M=testCase.Mechanical; S=testCase.Source;
            mapped=TestSegmentMechanicalMesh.interfaceCoordinates(M,S);
            expected=[M.p(M.interfaceRemapping.mechanicalEdgeNodes(:,1),:); ...
                M.p(M.interfaceRemapping.mechanicalEdgeNodes(:,2),:)];
            testCase.verifyTrue(M.interfaceRemapping.available);
            testCase.verifyEqual(mapped,expected,'AbsTol',1e-10);
            testCase.verifyFalse(testCase.Report.soilMeshChanged);
            testCase.verifyFalse(testCase.Report.thermalMeshChanged);
        end

        function eachSegmentCanMoveRigidlyRelativeToNeighbours(testCase)
            M=testCase.Mechanical;
            translation=[(1:8)'*1e-3,-((1:8)').^2*1e-3];
            moved=M.p+translation(M.nodeSegmentID,:);
            oldLengths=TestSegmentMechanicalMesh.quadEdgeLengths(M.p,M.q);
            newLengths=TestSegmentMechanicalMesh.quadEdgeLengths(moved,M.q);
            minus=reshape(M.jointNodePairs(:,1,:),[],1);
            plus=reshape(M.jointNodePairs(:,2,:),[],1);
            testCase.verifyEqual(newLengths,oldLengths,'AbsTol',1e-11);
            testCase.verifyTrue(all(hypot(moved(plus,1)-moved(minus,1), ...
                moved(plus,2)-moved(minus,2))>0));
        end

        function inputJointIDsSurvivePermutation(testCase)
            order=[8 1:7];
            M=build_segment_mechanical_mesh(testCase.Source,testCase.Angles(order));
            testCase.verifyEqual(M.p,testCase.Mechanical.p,'AbsTol',0);
            testCase.verifyEqual(M.q,testCase.Mechanical.q);
            testCase.verifyEqual(M.jointNodePairs(:,:,[2:8 1]), ...
                testCase.Mechanical.jointNodePairs);
        end

        function wraparoundJointAtZeroIsHandled(testCase)
            angles=(0:7)*2*pi/8;
            M=build_segment_mechanical_mesh(testCase.Source,angles);
            testCase.verifyEqual(numel(unique(M.connectedComponentID)),8);
            testCase.verifyEqual(M.jointSegmentPairs(1,:),[8 1]);
            testCase.verifyEqual(M.p(M.jointNodePairs(:,1,1),:), ...
                M.p(M.jointNodePairs(:,2,1),:),'AbsTol',1e-10);
        end

        function duplicatePhysicalJointAnglesAreRejected(testCase)
            angles=testCase.Angles;
            bad=[angles(1:7),angles(1)+2*pi];
            testCase.verifyError(@() build_segment_mechanical_mesh(testCase.Source,bad), ...
                'SegmentMechanicalMesh:JointAngles');
        end

        function invertBridgingCannotBeSilentlyEnabled(testCase)
            testCase.verifyError(@() build_segment_mechanical_mesh( ...
                testCase.Source,testCase.Angles,struct('includeInvert',true)), ...
                'SegmentMechanicalMesh:InvertUnsupported');
            testCase.verifyEqual(testCase.Report.excludedInvertQuadCount, ...
                size(testCase.Source.qInvertGlobal,1));
            testCase.verifyFalse(testCase.Mechanical.constitutiveLawsImplemented);
            testCase.verifyFalse(testCase.Mechanical.hydraulicFeedbackImplemented);
        end
    end

    methods (Static,Access=private)
        function components=graphComponents(q,N)
            a=[q(:,1);q(:,2);q(:,3);q(:,4)];
            b=[q(:,2);q(:,3);q(:,4);q(:,1)];
            components=conncomp(graph(a,b,[],N));
        end

        function lengths=quadEdgeLengths(p,q)
            x=reshape(p(q(:),1),size(q)); y=reshape(p(q(:),2),size(q));
            lengths=hypot(x-x(:,[2 3 4 1]),y-y(:,[2 3 4 1]));
        end

        function coordinates=interfaceCoordinates(M,S)
            I=M.interfaceRemapping; ne=size(I.mechanicalEdgeNodes,1);
            first=zeros(ne,2); last=first;
            for k=1:ne
                weights=squeeze(I.soilEndpointWeights(k,:,:));
                point=weights*S.pAll(I.soilGlobalNodes(k,:),:);
                first(k,:)=point(1,:); last(k,:)=point(2,:);
            end
            coordinates=[first;last];
        end

        function minimum=minimumQuadJacobian(p,q)
            minimum=inf;
            for k=1:size(q,1)
                for xi=[-1 1]/sqrt(3)
                    for eta=[-1 1]/sqrt(3)
                        [~,d]=TestSegmentMechanicalMesh.qshape(xi,eta);
                        minimum=min(minimum,det(d*p(q(k,:),:)));
                    end
                end
            end
        end

        function [child,source]=quadraturePullback(M,S,temperature)
            childTemperature=M.thermalInterpolation*temperature;
            nq=size(M.q,1); child=zeros(4*nq,1); source=child; row=0;
            for k=1:nq
                original=S.qTunnelGlobal(M.elementParentThermalQ4(k),:);
                for xi=[-1 1]/sqrt(3)
                    for eta=[-1 1]/sqrt(3)
                        row=row+1;
                        [N,~]=TestSegmentMechanicalMesh.qshape(xi,eta);
                        point=N*M.p(M.q(k,:),:);
                        shape=TestSegmentMechanicalMesh.inverseQuad(point,S.pAll(original,:));
                        child(row)=N*childTemperature(M.q(k,:));
                        source(row)=shape*temperature(original);
                    end
                end
            end
        end

        function shape=inverseQuad(point,xy)
            natural=[0;0];
            for iteration=1:12
                [shape,d]=TestSegmentMechanicalMesh.qshape(natural(1),natural(2));
                residual=(shape*xy-point)';
                natural=natural-(d*xy)'\residual;
            end
            [shape,~]=TestSegmentMechanicalMesh.qshape(natural(1),natural(2));
        end

        function [N,d]=qshape(xi,eta)
            N=.25*[(1-xi)*(1-eta),(1+xi)*(1-eta), ...
                (1+xi)*(1+eta),(1-xi)*(1+eta)];
            d=.25*[-(1-eta),1-eta,1+eta,-(1+eta); ...
                -(1-xi),-(1+xi),1+xi,1-xi];
        end
    end
end
