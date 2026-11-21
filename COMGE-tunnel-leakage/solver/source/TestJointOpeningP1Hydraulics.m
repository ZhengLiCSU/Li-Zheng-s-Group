classdef TestJointOpeningP1Hydraulics < matlab.unittest.TestCase
    % Endpoint continuity, exact physical geometry and coupled tangents.
    methods (Test)
        function uniformGapRecoversScalarThermalSeries(testCase)
            gap=repmat([.001;.002],1,3); weights=[.2 .3 .5];
            radii=[6 6;5.775 5.75;5.55 5.5];
            [B,U]=TestJointOpeningP1Hydraulics.fixture(gap,weights,radii);
            T=[-.01 .02 .7;-.008 .03 .5];
            options=struct('openingTransmissionFraction',[.001;.0015], ...
                'pathLength_m',[.45;.5],'freezeEnabled',true,'viscosityEnabled',true, ...
                'transitionWidth_C',.1,'sampleWeights',weights, ...
                'baseMobilityFactor',[.3;.4], ...
                'baseMobilityDerivative_dT',[.1 -.2 .3;-.2 .1 .4]);
            [expected,expectedG,expectedT]=joint_opening_hydraulics( ...
                [2e-10;3e-10],gap(:,1),T,options);
            [actual,~,actualT,D]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T,options);
            testCase.verifyEqual(actual,expected,'RelTol',2e-12,'AbsTol',1e-25);
            testCase.verifyEqual(sum(D.dC_dNodalGap,2),expectedG,'RelTol',2e-12,'AbsTol',1e-22);
            testCase.verifyEqual(actualT,expectedT,'RelTol',2e-12,'AbsTol',1e-25);
            testCase.verifyEqual(D.newslitGeometricVolume_m2, ...
                [.45;.5].*[.001;.0015].*gap(:,1),'AbsTol',1e-20);
        end

        function linearApertureHasAnalyticalCubicResistance(testCase)
            gap=[.001 .004]; fraction=.001;
            [B,U]=TestJointOpeningP1Hydraulics.fixture(gap,[.2 .3 .5]);
            options=struct('openingTransmissionFraction',fraction,'baseMobilityFactor',0);
            [actual,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[2 3 4],options);
            b=fraction*gap; L=B.radiusOuter_m-B.radiusInner_m;
            resistance=L*(b(1)+b(2))/(2*b(1)^2*b(2)^2);
            expected=1000*9.81/(12*water_viscosity_iapws(20)*resistance);
            testCase.verifyEqual(actual,expected,'RelTol',2e-13,'AbsTol',1e-26);
            testCase.verifyEqual(D.referenceOpenConductance,expected,'RelTol',2e-13,'AbsTol',1e-26);
            testCase.verifyEqual(D.openMobilityFactor,1,'AbsTol',1e-14);
        end

        function nonuniformRadialIntervalsRetainTheirActualLengths(testCase)
            gap=[.001 .004 .002]; fraction=.001; radii=[6;5.9;5.55];
            [B,U]=TestJointOpeningP1Hydraulics.fixture(gap,[.3 .7],radii);
            options=struct('openingTransmissionFraction',fraction,'baseMobilityFactor',0);
            actual=joint_opening_p1_hydraulics(2e-10,B,U,[3 4],options);
            b=fraction*gap;
            R=.1*(b(1)+b(2))/(2*b(1)^2*b(2)^2) ...
                +.35*(b(2)+b(3))/(2*b(2)^2*b(3)^2);
            expected=1000*9.81/(12*water_viscosity_iapws(20)*R);
            testCase.verifyEqual(actual,expected,'RelTol',3e-13,'AbsTol',1e-26);
        end

        function endpointClosingHasContinuousQuadraticLimit(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([1e-7 .002],[.5 .5]);
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',0);
            [coarse,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[4 4],options);
            Uf=U; Uf(1)=-9e-8;
            fine=joint_opening_p1_hydraulics(2e-10,B,Uf,[4 4],options);
            Uz=U; Uz(1)=-1e-7;
            [closed,closedU,closedT,Dclosed]=joint_opening_p1_hydraulics(2e-10,B,Uz,[4 4],options);
            testCase.verifyEqual(fine/coarse,.01,'RelTol',5e-5,'AbsTol',1e-8);
            testCase.verifyEqual(closed,0,'AbsTol',0);
            testCase.verifyEqual(full(closedU),zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(closedT,zeros(1,2),'AbsTol',0);
            testCase.verifyGreaterThan(D.dC_dNodalGap(1),0);
            testCase.verifyGreaterThan(Dclosed.newslitGeometricVolume_m2,0);
        end

        function midpointGapCannotCreateFiniteFlowAtClosedEndpoint(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([0 .002],[.5 .5]);
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',0);
            [actual,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[4 4],options);
            [~,~,~,midpoint]=joint_opening_profile_hydraulics( ...
                2e-10,[.0005 .0015],[4 4],options);
            Up=U; Up(1)=1e-9;
            nearOpen=joint_opening_p1_hydraulics(2e-10,B,Up,[4 4],options);
            testCase.verifyEqual(actual,0,'AbsTol',0);
            testCase.verifyFalse(D.wholePathOpen);
            testCase.verifyGreaterThan(midpoint.openConductance,0);
            testCase.verifyLessThan(nearOpen,1e-8*midpoint.openConductance);
        end

        function closedInteriorBlocksFlowAndRetainsExactInventory(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 -.00001 .003],[.5 .5]);
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',.4);
            [C,DU,DT,D]=joint_opening_p1_hydraulics(2e-10,B,U,[4 4],options);
            [expectedV,~,~,expectedDU]=joint_profile_geometric_volume(B,U,2e-10,.001);
            testCase.verifyEqual(C,.4*2e-10,'AbsTol',0);
            testCase.verifyEqual(full(DU),zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(DT,zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,expectedV,'AbsTol',1e-21);
            testCase.verifyEqual(D.dVolume_dU,expectedDU,'AbsTol',1e-16);
            testCase.verifyGreaterThan(expectedV,0);
        end

        function sealedJointHasNoNewPassageOrVolume(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .002 .003],[.5 .5]);
            options=struct('openingTransmissionFraction',1,'freezeEnabled',true,'transitionWidth_C',.1);
            [C,DU,DT,D]=joint_opening_p1_hydraulics(0,B,U,[-.01 .02],options);
            testCase.verifyEqual(C,0,'AbsTol',0);
            testCase.verifyEqual(full(DU),zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(DT,zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,0,'AbsTol',0);
            testCase.verifyTrue(D.sealedMask);
        end

        function iceBlockingPreservesGeometryAndIndependentBase(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .002 .003],[.5 .5]);
            options=struct('openingTransmissionFraction',.001,'freezeEnabled',true, ...
                'transitionWidth_C',.05,'baseMobilityFactor',.4);
            [C,DU,DT,D]=joint_opening_p1_hydraulics(2e-10,B,U,[4 -.1],options);
            testCase.verifyEqual(C,.4*2e-10,'AbsTol',0);
            testCase.verifyEqual(full(DU),zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(DT,zeros(1,2),'AbsTol',0);
            testCase.verifyGreaterThan(D.newslitGeometricVolume_m2,0);
            testCase.verifyGreaterThan(D.iceVolumeBySample_m2(2),0);
            testCase.verifyGreaterThan(D.referenceOpenConductance,0);
            testCase.verifyEqual(D.openMobilityFactor,0,'AbsTol',0);
        end

        function displacementAndTemperatureTangentsMatchDirectionalFD(testCase)
            [B,U,T,options]=TestJointOpeningP1Hydraulics.phaseFixture();
            directionU=[.13;-.07;.04;-.05;.08;.1];
            directionT=[.2 -.1 .3;.1 .15 -.2]; du=1e-8; dt=1e-6;
            [~,JU,JT]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T,options);
            pu=joint_opening_p1_hydraulics([2e-10;3e-10],B,U+du*directionU,T,options);
            mu=joint_opening_p1_hydraulics([2e-10;3e-10],B,U-du*directionU,T,options);
            pt=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T+dt*directionT,options);
            mt=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T-dt*directionT,options);
            testCase.verifyEqual(JU*directionU,(pu-mu)/(2*du),'RelTol',5e-7,'AbsTol',1e-17);
            testCase.verifyEqual(sum(JT.*directionT,2),(pt-mt)/(2*dt), ...
                'RelTol',5e-7,'AbsTol',1e-21);
        end

        function trueNodalVolumeAndPhaseInventoryTangentsMatchFD(testCase)
            [B,U,T,options]=TestJointOpeningP1Hydraulics.phaseFixture();
            vu=[.13;-.07;.04;-.05;.08;.1]; vt=[.2 -.1 .3;.1 .15 -.2]; du=1e-8; dt=1e-6;
            [~,~,~,D]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T,options);
            [~,~,~,pu]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U+du*vu,T,options);
            [~,~,~,mu]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U-du*vu,T,options);
            [~,~,~,pt]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T+dt*vt,options);
            [~,~,~,mt]=joint_opening_p1_hydraulics([2e-10;3e-10],B,U,T-dt*vt,options);
            testCase.verifyEqual(D.dVolume_dU*vu, ...
                (pu.newslitGeometricVolume_m2-mu.newslitGeometricVolume_m2)/(2*du), ...
                'RelTol',1e-8,'AbsTol',1e-13);
            testCase.verifyEqual(D.dLiquidVolume_dU*vu,(pu.liquidVolume_m2-mu.liquidVolume_m2)/(2*du), ...
                'RelTol',1e-8,'AbsTol',1e-13);
            testCase.verifyEqual(D.dLiquidVolume_dTemperature_m2_K.*vt, ...
                (pt.liquidVolumeBySample_m2-mt.liquidVolumeBySample_m2)/(2*dt), ...
                'RelTol',1e-7,'AbsTol',1e-15);
            testCase.verifyEqual(D.liquidVolumeBySample_m2+D.iceVolumeBySample_m2, ...
                D.sampleVolume_m2,'AbsTol',1e-21);
        end

        function falsePathMaskPreservesInventoryAndBaseDerivative(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .002 .003],[.2 .3 .5]);
            options=struct('openingTransmissionFraction',.001,'openPathMask',false, ...
                'baseMobilityFactor',.4,'baseMobilityDerivative_dT',[.1 -.2 .3]);
            [C,DU,DT,D]=joint_opening_p1_hydraulics(2e-10,B,U,[2 3 4],options);
            testCase.verifyEqual(C,.4*2e-10,'AbsTol',0);
            testCase.verifyEqual(full(DU),zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(DT,2e-10*[.1 -.2 .3],'AbsTol',1e-25);
            testCase.verifyGreaterThan(D.newslitGeometricVolume_m2,0);
            testCase.verifyGreaterThan(D.referenceOpenConductance,0);
        end

        function reversedNodalOrderRetainsConductanceAndDerivativeIdentity(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .003 .002],[.2 .3 .5]);
            [Br,Ur]=TestJointOpeningP1Hydraulics.fixture([.002 .003 .001], ...
                [.2 .3 .5],flipud(B.nodalRadius_m));
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',0);
            [C,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[2 3 4],options);
            [Cr,~,~,Dr]=joint_opening_p1_hydraulics(2e-10,Br,Ur,[2 3 4],options);
            testCase.verifyEqual(Cr,C,'RelTol',1e-13,'AbsTol',1e-26);
            testCase.verifyEqual(Dr.dC_dNodalGap,fliplr(D.dC_dNodalGap),'RelTol',1e-13,'AbsTol',1e-22);
        end

        function refiningTheSameLinearProfileLeavesFlowAndVolumeUnchanged(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .004],[.2 .3 .5]);
            [Bf,Uf]=TestJointOpeningP1Hydraulics.fixture([.001 .0025 .004],[.2 .3 .5]);
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',0);
            [C,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[2 3 4],options);
            [Cf,~,~,Df]=joint_opening_p1_hydraulics(2e-10,Bf,Uf,[2 3 4],options);
            testCase.verifyEqual(Cf,C,'RelTol',3e-13,'AbsTol',1e-26);
            testCase.verifyEqual(Df.newslitGeometricVolume_m2,D.newslitGeometricVolume_m2,'AbsTol',1e-21);
            testCase.verifyEqual(sum(Df.dC_dNodalGap,2),sum(D.dC_dNodalGap,2), ...
                'RelTol',3e-13,'AbsTol',1e-22);
        end

        function unhitNarrowPositiveEndpointStillHasExactVolume(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.00002 -.00018],[.5 .5]);
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',0);
            [C,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[2 3],options);
            [V,~,~,DU]=joint_profile_geometric_volume(B,U,2e-10,.001);
            testCase.verifyEqual(C,0,'AbsTol',0);
            testCase.verifyGreaterThan(V,0);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,V,'AbsTol',1e-23);
            testCase.verifyEqual(D.dVolume_dU,DU,'AbsTol',1e-16);
        end

        function tinyPositiveFlowDoesNotNeedAnApertureFloor(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([1e-85 2e-85],[.5 .5]);
            options=struct('openingTransmissionFraction',1,'baseMobilityFactor',0);
            [C,DU,~,D]=joint_opening_p1_hydraulics(2e-10,B,U,[2 3],options);
            testCase.verifyGreaterThan(C,0);
            testCase.verifyLessThan(C,1e-240);
            testCase.verifyGreaterThan(full(DU),zeros(1,2));
            testCase.verifyTrue(all(isfinite(D.dC_dNodalGap)));
        end

        function apertureDependentFreezingRequiresExplicitSpatialTreatment(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .002],[.5 .5]);
            options=struct('openingTransmissionFraction',.001,'wallEnergyDifference_J_m2',.015);
            testCase.verifyError(@()joint_opening_p1_hydraulics(2e-10,B,U,[2 3],options), ...
                'JointOpeningP1:ApertureDependentFreezing');
        end

        function anOlderCacheCanBeUsedAndDefaultMappingRemainsDisabled(testCase)
            [B,U]=TestJointOpeningP1Hydraulics.fixture([.001 .002],[.5 .5]);
            old=rmfield(B,'resistanceIntegration');
            [C,DU,DT,D]=joint_opening_p1_hydraulics(2e-10,old,U,[],struct());
            testCase.verifyEqual(C,2e-10,'AbsTol',0);
            testCase.verifyEqual(full(DU),zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(DT,zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,0,'AbsTol',0);
            testCase.verifyFalse(D.thermalProfileProvided);
        end
    end

    methods (Static,Access=private)
        function [B,U]=fixture(gap,weights,radius)
            [nj,nr]=size(gap); nq=numel(weights);
            if nargin<3, radius=repmat(linspace(6,5.55,nr)',1,nj); end
            B=struct('nj',nj,'nr',nr,'nq',nq,'nodalRadius_m',radius, ...
                'nodalGap0',reshape(gap',[],1),'nodalGapU',speye(nj*nr), ...
                'radiusOuter_m',max(radius,[],1)','radiusInner_m',min(radius,[],1)', ...
                'sampleWeights',weights);
            B.resistanceIntegration=joint_opening_p1_hydraulics('build',B);
            U=zeros(nj*nr,1);
        end

        function [B,U,T,options]=phaseFixture()
            gap=[.001 .0018 .003;.0012 .0022 .0035];
            [B,U]=TestJointOpeningP1Hydraulics.fixture(gap,[.2 .3 .5]);
            T=[.004 .015 .7;.012 .018 .5];
            options=struct('openingTransmissionFraction',[.001;.0015], ...
                'freezeEnabled',true,'viscosityEnabled',true,'transitionWidth_C',.05, ...
                'baseMobilityFactor',[.2;.4]);
        end
    end
end
