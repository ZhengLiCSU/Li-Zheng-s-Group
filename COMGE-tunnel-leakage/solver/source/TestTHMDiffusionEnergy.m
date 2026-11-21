classdef TestTHMDiffusionEnergy < matlab.unittest.TestCase
    %TESTTHMDIFFUSIONENERGY Physical P1/CVFE diffusion on obtuse triangles.
    % Pure mathematical/public-kernel fixtures, not a tunnel case or field
    % validation. Negative aggregate cotangent coefficients are retained.

    methods (Test)
        function negativeCotangentStillHasNonnegativeElementEnergy(testCase)
            r=TestTHMDiffusionEnergy.fixture(false,false);
            testCase.verifyGreaterThan(r.negativeEdgeCount,0);
            testCase.verifyGreaterThan(r.K(1,2),0);
            testCase.verifyGreaterThan(r.elementEnergy,0);
            testCase.verifyEqual(r.matrixEnergy,r.elementEnergy,'AbsTol',1e-19);
            testCase.verifyGreaterThanOrEqual(min(eig(full(r.K))),-1e-18);
            testCase.verifyGreaterThan(min(eig(full(r.K(2:end,2:end)))),0);
            testCase.verifyEqual(norm(r.K-r.K','inf'),0,'AbsTol',1e-20);
        end

        function physicalSubfaceOperatorEqualsGalerkinDiffusion(testCase)
            r=TestTHMDiffusionEnergy.fixture(false,false);
            testCase.verifyEqual(full(r.cv.inc*r.DH),full(r.K),'AbsTol',1e-20);
            testCase.verifyEqual(r.cv.inc*r.F,r.K*r.H,'AbsTol',1e-20);
            testCase.verifyEqual(r.directDarcySubface,r.F,'AbsTol',1e-20);
        end

        function constantPotentialAndGravityBalancedPatchAreExact(testCase)
            r=TestTHMDiffusionEnergy.fixture(false,false);
            [flow,derivative]=soil_cvfe_flux(r.cv.z+2,r.Ke,r.cv,'element');
            testCase.verifyEqual(flow,zeros(size(flow)),'AbsTol',1e-20);
            testCase.verifyEqual(r.K*ones(r.cv.ns,1),zeros(r.cv.ns,1),'AbsTol',1e-20);
            testCase.verifyEqual(r.cv.inc*derivative*ones(r.cv.ns,1),zeros(r.cv.ns,1),'AbsTol',1e-20);
        end

        function unitKPositiveAggregationDoesNotProveHeterogeneousSign(testCase)
            r=TestTHMDiffusionEnergy.fixture(true,false);
            testCase.verifyGreaterThanOrEqual(min(r.unitEdgeWeights),0);
            testCase.verifyGreaterThan(r.negativeEdgeCount,0);
            testCase.verifyLessThan(min(r.singleElementEdgeWeights),0);
            testCase.verifyGreaterThan(r.elementEnergy,0);
            testCase.verifyEqual(r.matrixEnergy,r.elementEnergy,'AbsTol',1e-18);
            testCase.verifyEqual(full(r.cv.inc*r.DH),full(r.K),'AbsTol',1e-20);
        end

        function clockwiseTriangleUsesTheSamePhysicalDiffusion(testCase)
            forward=TestTHMDiffusionEnergy.fixture(false,false);
            reverse=TestTHMDiffusionEnergy.fixture(false,true);
            testCase.verifyEqual(full(reverse.K),full(forward.K),'AbsTol',1e-20);
            testCase.verifyEqual(full(reverse.cv.inc*reverse.DH),full(forward.K),'AbsTol',1e-20);
            testCase.verifyEqual(reverse.cv.inc*reverse.F,forward.cv.inc*forward.F,'AbsTol',1e-20);
            testCase.verifyEqual(reverse.cv.volume,forward.cv.volume,'AbsTol',1e-15);
        end

        function heatDonorFollowsPhysicalFaceRatherThanSignedCotangent(testCase)
            r=TestTHMDiffusionEnergy.fixture(false,false);
            [heat,~,~,diagnostic]=soil_cvfe_heat([0;10;20],r.F,r.cv,1000,3.34e5,4180);
            face=find(r.cv.from==1&r.cv.to==2);
            testCase.verifyLessThan(r.F(face),0);
            testCase.verifyGreaterThan(r.aggregateEdge12Flow,0);
            testCase.verifyEqual(diagnostic.upstreamGlobalNode(face),2);
            testCase.verifyEqual(sum(heat),0,'AbsTol',1e-10);
        end
    end

    methods (Static,Access=private)
        function r=fixture(heterogeneous,reverse)
            xy=[0 0;2 0;.3 .2];tri=[1 2 3];Ke=1e-7;H=[0;1;1];
            if heterogeneous
                xy=[xy;1 -3];tri=[1 2 3;2 1 4];Ke=[1e-6;1e-7];H=[0;1;1;-.4];
            end
            if reverse,tri=tri(:,[1 3 2]);end
            nodes=(1:size(xy,1))';cv=build_soil_cvfe_flux(xy,tri,nodes);
            [flow,DH]=soil_cvfe_flux(H+cv.z,Ke,cv,'element');
            ne=size(tri,1);ns=size(xy,1);K=sparse(ns,ns);single=zeros(ne,3);
            gradients=zeros(ne,2);
            for e=1:ne
                points=xy(tri(e,:),:);twiceArea=det([points(2,:)-points(1,:);points(3,:)-points(1,:)]);
                grad=[points(2,2)-points(3,2),points(3,1)-points(2,1); ...
                    points(3,2)-points(1,2),points(1,1)-points(3,1); ...
                    points(1,2)-points(2,2),points(2,1)-points(1,1)]/twiceArea;
                stiffness=abs(twiceArea)/2*(grad*grad');
                K(tri(e,:),tri(e,:))=K(tri(e,:),tri(e,:))+Ke(e)*stiffness;
                single(e,:)=-[stiffness(1,2),stiffness(2,3),stiffness(3,1)];
                gradients(e,:)=H(tri(e,:))'*grad;
            end
            edgePairs=[tri(:,[1 2]);tri(:,[2 3]);tri(:,[3 1])];
            [edges,~,map]=unique(sort(edgePairs,2),'rows');
            unit=accumarray(map,single(:),[size(edges,1) 1]);
            weighted=accumarray(map,single(:).*repmat(Ke,3,1),[size(edges,1) 1]);
            edge12=find(edges(:,1)==1&edges(:,2)==2);
            direct=-Ke(cv.elementIndex).*sum(gradients(cv.elementIndex,:).*cv.normalLength,2);
            r=struct('K',K,'Ke',Ke,'H',H,'cv',cv,'F',flow,'DH',DH, ...
                'elementEnergy',sum(cv.area.*Ke.*sum(gradients.^2,2)), ...
                'matrixEnergy',H'*(K*H),'directDarcySubface',direct, ...
                'negativeEdgeCount',sum(weighted<0),'singleElementEdgeWeights',single(:), ...
                'unitEdgeWeights',unit,'aggregateEdge12Flow',weighted(edge12)*(H(1)-H(2)));
        end
    end
end
