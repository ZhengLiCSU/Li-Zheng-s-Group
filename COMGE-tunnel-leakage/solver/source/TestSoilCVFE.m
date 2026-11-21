classdef TestSoilCVFE < matlab.unittest.TestCase
    %TESTSOILCVFE Physical Darcy subfaces, P1 equivalence and heat tangents.
    % Independent small meshes need no tunnel files or nonlinear solver.
    % The stiffness oracle integrates grad(Ni).grad(Nj) directly; it does
    % not construct a median dual or use the kernel incidence matrix.

    properties (TestParameter)
        orientation=struct('counterclockwise',[1 2 3;2 4 3], ...
            'clockwise',[3 2 1;3 4 2]);
    end

    methods (Test)
        function obtuseTriangleUsesPhysicalNormalsDespiteNegativeCotangent(testCase)
            xy=[0 0;2 0;.5 .2]; tri=[1 2 3];
            G=build_soil_cvfe_flux(xy,tri,(1:3)');
            K=2; h=xy(:,1)+xy(:,2);
            expectedNormal=[1/15 1/6;-1/30 5/12;-1/30 -7/12];
            A=TestSoilCVFE.directStiffness(xy,tri,K);
            legacyPairFlow=-A(1,2)*(xy(1,1)-xy(2,1));
            [F,~,~,d]=soil_cvfe_flux(h,K,G,'element');
            testCase.verifyGreaterThan(A(1,2),0);
            testCase.verifyGreaterThan(legacyPairFlow,0);
            testCase.verifyLessThan(F(1),0);
            testCase.verifyEqual(G.normalLength,expectedNormal,'AbsTol',1e-14);
            testCase.verifyEqual(d.darcyElement,[-K 0],'AbsTol',1e-13);
            testCase.verifyEqual(F,expectedNormal*[-K;0],'AbsTol',1e-13);
            testCase.verifyEqual(G.inc*F,A*xy(:,1),'AbsTol',1e-13);
        end

        function affineHeadRecoversDarcyVectorAndP1ForEitherOrientation(testCase,orientation)
            xy=[0 0;1.3 .2;-.2 .9;1.1 1.3];
            G=build_soil_cvfe_flux(xy,orientation,(1:4)');
            gradient=[.6;-1.7]; K=3.2;
            H=xy*gradient+.7; h=H+xy(:,2);
            [F,DH,~,d]=soil_cvfe_flux(h,K*ones(4,1),G);
            A=TestSoilCVFE.directStiffness(xy,orientation,K*ones(2,1));
            expectedQ=repmat((-K*gradient)',2,1);
            expectedF=sum(G.normalLength.*expectedQ(G.elementIndex,:),2);
            testCase.verifyEqual(d.darcyElement,expectedQ,'AbsTol',2e-14);
            testCase.verifyEqual(F,expectedF,'AbsTol',3e-14);
            testCase.verifyEqual(full(G.inc*DH),A,'AbsTol',3e-14);
            testCase.verifyEqual(G.inc*F,A*H,'AbsTol',3e-14);
        end

        function varyingNodalKMatchesIndependentElementMeanStiffness(testCase)
            [xy,tri,soilNodes,G]=TestSoilCVFE.meshFixture();
            h=[1.3;-.4;2.1;.7]; K=[.2;1.7;.8;2.3];
            Ke=mean(K(G.localTri),2);
            [F,DH,DK]=soil_cvfe_flux(h,K,G);
            [Fe,DHe,DKe]=soil_cvfe_flux(h,Ke,G,'element');
            A=TestSoilCVFE.directStiffness(xy(soilNodes,:),G.localTri,Ke);
            testCase.verifyEqual(G.triGlobal,tri);
            testCase.verifyEqual(F,Fe,'AbsTol',3e-14);
            testCase.verifyEqual(full(DH),full(DHe),'AbsTol',3e-14);
            testCase.verifyEqual(full(DK),full(DKe*G.meanK),'AbsTol',3e-14);
            testCase.verifyEqual(full(G.inc*DH),A,'AbsTol',3e-14);
            testCase.verifyEqual(G.inc*F,A*(h-G.z),'AbsTol',3e-14);
            testCase.verifyEqual(sum(G.volume),sum(G.area),'AbsTol',1e-14);
        end

        function duplicatePairSubfacesRemainDistinctBeforeUpwinding(testCase)
            xy=[0 0;1 0;0 1;1.7 1.2]; tri=[1 2 3;2 4 3];
            G=build_soil_cvfe_flux(xy,tri,(1:4)');
            h=[0;1;0;10]+xy(:,2); K=[1;2];
            F=soil_cvfe_flux(h,K,G,'element');
            T=[1;2;5;8]; [~,~,~,d]=soil_cvfe_heat(T,F,G,1,0,1);
            testCase.verifySize(F,[6 1]);
            testCase.verifyEqual(G.from([3 6]),[2;3]);
            testCase.verifyEqual(G.to([3 6]),[3;2]);
            testCase.verifyGreaterThan(F(3),0);
            testCase.verifyGreaterThan(F(6),0);
            testCase.verifyEqual(d.upstreamGlobalNode([3 6]),[2;3]);
        end

        function hydrostaticPressureHasNoFlux(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            [F,DH]=soil_cvfe_flux(G.z+4,ones(G.ns,1),G);
            testCase.verifyEqual(F,zeros(G.nf,1),'AbsTol',2e-14);
            testCase.verifyEqual(DH*ones(G.ns,1),zeros(G.nf,1),'AbsTol',1e-14);
        end

        function zeroMobilityStopsWaterAndHeatWithoutFloor(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            [F,DH]=soil_cvfe_flux([1;2;3;4],zeros(G.ne,1),G,'element');
            [r,DT]=soil_cvfe_heat((1:G.N)',F,G,1000,334000,4180);
            testCase.verifyEqual(F,zeros(G.nf,1),'AbsTol',0);
            testCase.verifyEqual(full(DH),zeros(G.nf,G.ns),'AbsTol',0);
            testCase.verifyEqual(r,zeros(G.N,1),'AbsTol',0);
            testCase.verifyEqual(full(DT),zeros(G.N),'AbsTol',0);
        end

        function fluxHeadAndNodalConductivityTangentsMatchFiniteDifference(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            h=[1.3;-.4;2.1;.7]; K=[.2;1.7;.8;2.3];
            vh=[.3;-.1;.5;-.2]; vK=[-.02;.03;-.04;.05]; e=1e-6;
            [~,DH,DK]=soil_cvfe_flux(h,K,G);
            ph=soil_cvfe_flux(h+e*vh,K,G); mh=soil_cvfe_flux(h-e*vh,K,G);
            pk=soil_cvfe_flux(h,K+e*vK,G); mk=soil_cvfe_flux(h,K-e*vK,G);
            testCase.verifyEqual(DH*vh,(ph-mh)/(2*e),'AbsTol',1e-9);
            testCase.verifyEqual(DK*vK,(pk-mk)/(2*e),'AbsTol',1e-9);
            testCase.verifyTrue(issparse(DH));
            testCase.verifyTrue(issparse(DK));
        end

        function elementConductivityTangentMatchesFiniteDifference(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            h=[1.3;-.4;2.1;.7]; Ke=[.4;1.9]; v=[.07;-.13]; e=1e-6;
            [~,~,DK]=soil_cvfe_flux(h,Ke,G,'element');
            p=soil_cvfe_flux(h,Ke+e*v,G,'element');
            m=soil_cvfe_flux(h,Ke-e*v,G,'element');
            testCase.verifyEqual(DK*v,(p-m)/(2*e),'AbsTol',1e-9);
            testCase.verifySize(DK,[G.nf G.ne]);
        end

        function globalHeatConservationIncludesTangentColumnSums(testCase)
            [~,~,soilNodes,G]=TestSoilCVFE.meshFixture();
            T=[-2;1;5;.7;20;-30]; F=[.2;-.5;.7;-.3;-.2;.4]*1e-7;
            [r,DT,DF]=soil_cvfe_heat(T,F,G,1000,334000,4180);
            testCase.verifyEqual(sum(r),0,'AbsTol',2e-13);
            testCase.verifyEqual(full(sum(DT,1)),zeros(1,G.N),'AbsTol',1e-14);
            testCase.verifyEqual(full(sum(DF,1)),zeros(1,G.nf),'AbsTol',1e-7);
            testCase.verifyEqual(r(setdiff((1:G.N)',soilNodes)),zeros(2,1),'AbsTol',0);
            testCase.verifyEqual(full(DT(:,setdiff((1:G.N)',soilNodes))), ...
                zeros(G.N,2),'AbsTol',0);
        end

        function upstreamChoiceFollowsPhysicalFaceFlow(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            T=(1:G.N)'; F=[2;-3;4;-5;6;-7]*1e-8;
            [r,~,~,d]=soil_cvfe_heat(T,F,G,2,3,5);
            expectedDonor=[G.fromGlobal(1);G.toGlobal(2);G.fromGlobal(3); ...
                G.toGlobal(4);G.fromGlobal(5);G.toGlobal(6)];
            expectedFlow=2*(3+5*T(expectedDonor)).*F;
            testCase.verifyEqual(d.upstreamGlobalNode,expectedDonor);
            testCase.verifyEqual(d.advectiveEnthalpyFlux,expectedFlow,'AbsTol',1e-20);
            testCase.verifyEqual(r,G.soilSelect'*G.inc*expectedFlow,'AbsTol',1e-20);
        end

        function heatTemperatureAndFlowTangentsMatchFiniteDifference(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            T=[-2;1;5;.7;20;-30]; F=[.2;-.5;.7;-.3;-.2;.4]*1e-7;
            vT=[.2;-.1;.4;.3;-.2;.1]; vF=[.1;-.2;.3;.2;-.1;.4]*1e-7;
            e=1e-5; [~,DT,DF]=soil_cvfe_heat(T,F,G,1000,334000,4180);
            pT=soil_cvfe_heat(T+e*vT,F,G,1000,334000,4180);
            mT=soil_cvfe_heat(T-e*vT,F,G,1000,334000,4180);
            pF=soil_cvfe_heat(T,F+e*vF,G,1000,334000,4180);
            mF=soil_cvfe_heat(T,F-e*vF,G,1000,334000,4180);
            testCase.verifyEqual(DT*vT,(pT-mT)/(2*e),'AbsTol',3e-9);
            testCase.verifyEqual(DF*vF,(pF-mF)/(2*e),'AbsTol',3e-9);
        end

        function nonlinearElementMobilityChainTangentMatchesFiniteDifference(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            h=[1.3;-.4;2.1;.7]; T=[-2;1;5;.7;20;-30];
            vH=[.3;-.1;.5;-.2]; vT=[.2;-.1;.4;.3;-.2;.1];
            [~,JT,JH]=TestSoilCVFE.compositeHeat(h,T,G);
            e=1e-5;
            pT=TestSoilCVFE.compositeHeat(h,T+e*vT,G);
            mT=TestSoilCVFE.compositeHeat(h,T-e*vT,G);
            pH=TestSoilCVFE.compositeHeat(h+e*vH,T,G);
            mH=TestSoilCVFE.compositeHeat(h-e*vH,T,G);
            testCase.verifyEqual(JT*vT,(pT-mT)/(2*e),'AbsTol',3e-9);
            testCase.verifyEqual(JH*vH,(pH-mH)/(2*e),'AbsTol',3e-9);
        end

        function invalidTriangleAndConductivityAreRejected(testCase)
            [~,~,~,G]=TestSoilCVFE.meshFixture();
            testCase.verifyError(@()build_soil_cvfe_flux([0 0;1 0;2 0],[1 2 3]), ...
                'build_soil_cvfe_flux:DegenerateTriangle');
            testCase.verifyError(@()soil_cvfe_flux(ones(4,1),[1;2;3],G,'element'), ...
                'soil_cvfe_flux:ConductivitySize');
        end
    end

    methods (Static,Access=private)
        function [xy,tri,soilNodes,G]=meshFixture()
            xy=[0 0;1.3 .2;-.2 .9;1.1 1.3;3 2;4 2];
            tri=[1 2 3;2 4 3]; soilNodes=[3;1;4;2];
            G=build_soil_cvfe_flux(xy,tri,soilNodes);
        end

        function A=directStiffness(xy,tri,Ke)
            A=zeros(size(xy,1));
            for element=1:size(tri,1)
                ids=tri(element,:); vertices=xy(ids,:);
                coefficients=[ones(3,1),vertices]\eye(3);
                gradients=coefficients(2:3,:);
                area=abs(det([vertices(2,:)-vertices(1,:); ...
                    vertices(3,:)-vertices(1,:)]))/2;
                A(ids,ids)=A(ids,ids)+Ke(element)*area*(gradients'*gradients);
            end
        end

        function [r,JT,JH]=compositeHeat(h,T,G)
            % Independent smooth law only for differentiating composition.
            localT=G.soilSelect*T;
            Ke=1e-8*exp(.01*(G.meanK*localT)+.02*(G.meanK*h));
            dKeT=spdiags(.01*Ke,0,G.ne,G.ne)*G.meanK*G.soilSelect;
            dKeh=spdiags(.02*Ke,0,G.ne,G.ne)*G.meanK;
            [F,DH,DK]=soil_cvfe_flux(h,Ke,G,'element');
            [r,DT,DF]=soil_cvfe_heat(T,F,G,1000,334000,4180);
            JT=DT+DF*DK*dKeT;
            JH=DF*(DH+DK*dKeh);
        end
    end
end
