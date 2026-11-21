classdef TestJointOpeningProfileHydraulics < matlab.unittest.TestCase
    %TESTJOINTOPENINGPROFILEHYDRAULICS Variable-aperture series/volume checks.
    % Constitutive benchmarks require no mesh, THM core or site data.

    methods (Test)
        function uniformApertureRecoversScalarKernelWithThermalProfile(testCase)
            cref=[2e-10;3e-10]; gap=[.001;.002];
            T=[-.006 .02 .7;-.008 .03 .5];
            options=struct('openingTransmissionFraction',[.001;.0015], ...
                'pathLength_m',[.45;.5],'freezeEnabled',true,'viscosityEnabled',true, ...
                'wallEnergyDifference_J_m2',.015,'transitionWidth_C',.1, ...
                'sampleWeights',[.2 .3 .5],'baseMobilityFactor',[.3;.4], ...
                'baseMobilityDerivative_dT',[.1 -.2 .3;-.2 .1 .4]);
            [C0,Dg0,Dt0,D0]=joint_opening_hydraulics(cref,gap,T,options);
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(cref,repmat(gap,1,3),T,options);
            testCase.verifyEqual(C,C0,'RelTol',2e-13,'AbsTol',1e-25);
            testCase.verifyEqual(sum(Dg,2),Dg0,'RelTol',2e-12,'AbsTol',1e-20);
            testCase.verifyEqual(Dt,Dt0,'RelTol',2e-12,'AbsTol',1e-24);
            testCase.verifyEqual(D.openConductance,D0.openConductance,'RelTol',2e-13);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,D0.newslitGeometricVolume_m2,'AbsTol',1e-21);
            testCase.verifyEqual(sum(D.dVolume_dGap_m,2),D0.dVolume_dGap_m,'AbsTol',1e-18);
            testCase.verifyEqual(D.liquidVolumeBySample_m2,D0.liquidVolumeBySample_m2,'AbsTol',1e-21);
        end

        function variableApertureMatchesAnalyticalPiecewiseSlitResistance(testCase)
            gap=[.001 .002 .004]; weights=[.2 .3 .5]; fraction=.001;
            options=struct('openingTransmissionFraction',fraction,'sampleWeights',weights);
            coefficient=1000*9.81/(12*TestJointOpeningProfileHydraulics.viscosity(20)*.45);
            b=fraction*gap; expected=coefficient/sum(weights./b.^3);
            [C,Dg,~,D]=joint_opening_profile_hydraulics(2e-10,gap,[2 3 4],options);
            expectedSlope=expected^2*weights.*(3./(coefficient*fraction^3*gap.^4));
            testCase.verifyEqual(D.openConductance,expected,'RelTol',1e-13);
            testCase.verifyEqual(C,2e-10+expected,'RelTol',1e-13);
            testCase.verifyEqual(Dg,expectedSlope,'RelTol',1e-13);
            testCase.verifyGreaterThan(D.openConductance,min(D.localConductance));
            testCase.verifyLessThan(D.openConductance,max(D.localConductance));
        end

        function anyClosedPointBlocksOnlyTheNewSeriesPassage(testCase)
            options=struct('openingTransmissionFraction',.001,'baseMobilityFactor',.4);
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(2e-10,[.001 0 .003],[2 3 4],options);
            testCase.verifyEqual(C,.4*2e-10,'AbsTol',0);
            testCase.verifyEqual(D.openConductance,0,'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(Dt,zeros(1,3),'AbsTol',0);
            testCase.verifyGreaterThan(D.newslitGeometricVolume_m2,0);
        end

        function localFreezingCanBlockAnOtherwiseOpenProfile(testCase)
            options=struct('openingTransmissionFraction',.001,'freezeEnabled',true);
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(2e-10,[.001 .002 .003],[2 -.1 4],options);
            testCase.verifyEqual(C,2e-10,'AbsTol',0);
            testCase.verifyEqual(D.openConductance,0,'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(Dt,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(D.liquidFractionBySample,[1 0 1],'AbsTol',0);
            testCase.verifyGreaterThan(D.iceVolumeBySample_m2(2),0);
        end

        function sealedJointCannotCreatePassageOrVolume(testCase)
            options=struct('openingTransmissionFraction',1,'freezeEnabled',true, ...
                'transitionWidth_C',.1,'wallEnergyDifference_J_m2',.015);
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(0,[.001 .002 .003],[-.01 .02 .7],options);
            testCase.verifyEqual(C,0,'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(Dt,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,0,'AbsTol',0);
            testCase.verifyEqual(D.sampleVolume_m2,zeros(1,3),'AbsTol',0);
            testCase.verifyTrue(D.sealedMask);
        end

        function independentlyVaryingGapAndTemperatureTangentsMatchFD(testCase)
            [cref,gap,T,options]=TestJointOpeningProfileHydraulics.phaseFixture();
            vg=[.13 -.07 .04;-.05 .08 .1]; vt=[.2 -.1 .3;.1 .15 -.2];
            eg=1e-8; et=1e-6;
            [~,Dg,Dt]=joint_opening_profile_hydraulics(cref,gap,T,options);
            pg=joint_opening_profile_hydraulics(cref,gap+eg*vg,T,options);
            mg=joint_opening_profile_hydraulics(cref,gap-eg*vg,T,options);
            pt=joint_opening_profile_hydraulics(cref,gap,T+et*vt,options);
            mt=joint_opening_profile_hydraulics(cref,gap,T-et*vt,options);
            testCase.verifyEqual(sum(Dg.*vg,2),(pg-mg)/(2*eg),'RelTol',5e-7,'AbsTol',1e-15);
            testCase.verifyEqual(sum(Dt.*vt,2),(pt-mt)/(2*et),'RelTol',5e-7,'AbsTol',1e-18);
        end

        function weightedGeometricAndPhaseInventoriesDifferentiateCorrectly(testCase)
            [cref,gap,T,options]=TestJointOpeningProfileHydraulics.phaseFixture();
            vg=[.13 -.07 .04;-.05 .08 .1]; vt=[.2 -.1 .3;.1 .15 -.2];
            eg=1e-8; et=1e-6;
            [~,~,~,D]=joint_opening_profile_hydraulics(cref,gap,T,options);
            [~,~,~,pg]=joint_opening_profile_hydraulics(cref,gap+eg*vg,T,options);
            [~,~,~,mg]=joint_opening_profile_hydraulics(cref,gap-eg*vg,T,options);
            [~,~,~,pt]=joint_opening_profile_hydraulics(cref,gap,T+et*vt,options);
            [~,~,~,mt]=joint_opening_profile_hydraulics(cref,gap,T-et*vt,options);
            expectedV=options.pathLength_m.*sum(gap.*options.openingTransmissionFraction.*options.sampleWeights,2);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,expectedV,'AbsTol',1e-21);
            testCase.verifyEqual(D.liquidVolumeBySample_m2+D.iceVolumeBySample_m2, ...
                D.sampleVolume_m2,'AbsTol',1e-21);
            testCase.verifyEqual(sum(D.dVolume_dGap_m.*vg,2), ...
                (pg.newslitGeometricVolume_m2-mg.newslitGeometricVolume_m2)/(2*eg), ...
                'RelTol',1e-8,'AbsTol',1e-12);
            testCase.verifyEqual(D.dLiquidVolume_dGap_m.*vg, ...
                (pg.liquidVolumeBySample_m2-mg.liquidVolumeBySample_m2)/(2*eg), ...
                'RelTol',5e-7,'AbsTol',1e-12);
            testCase.verifyEqual(D.dLiquidVolume_dTemperature_m2_K.*vt, ...
                (pt.liquidVolumeBySample_m2-mt.liquidVolumeBySample_m2)/(2*et), ...
                'RelTol',5e-7,'AbsTol',1e-15);
        end

        function everyPointHasItsOwnGTFreezingTemperature(testCase)
            gap=[.001 .002 .004]; fraction=.001; gamma=.015;
            options=struct('openingTransmissionFraction',fraction,'freezeEnabled',true, ...
                'wallEnergyDifference_J_m2',gamma);
            [~,~,~,D]=joint_opening_profile_hydraulics(2e-10,gap,[1 2 3],options);
            expected=-2*273.15*gamma./(917*3.35e5*fraction*gap);
            testCase.verifyEqual(D.freezeTemperature_C,expected,'RelTol',1e-14);
            testCase.verifyGreaterThan(diff(D.freezeTemperature_C),zeros(1,2));
        end

        function singleSampleIsExactlyScalarKernel(testCase)
            cref=[2e-10;3e-10]; gap=[.001;.002]; T=[2;4];
            options=struct('openingTransmissionFraction',.001,'viscosityEnabled',true);
            [C0,Dg0,Dt0,D0]=joint_opening_hydraulics(cref,gap,T,options);
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(cref,gap,T,options);
            testCase.verifyEqual(C,C0,'RelTol',1e-14);
            testCase.verifyEqual(Dg,Dg0,'RelTol',1e-14);
            testCase.verifyEqual(Dt,Dt0,'RelTol',1e-14);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,D0.newslitGeometricVolume_m2,'AbsTol',0);
        end

        function invalidTemperatureShapeAndZeroWeightAreRejected(testCase)
            testCase.verifyError(@()joint_opening_profile_hydraulics(2e-10, ...
                [.001 .002],[2;3],struct()),'JointOpeningProfile:Temperature');
            testCase.verifyError(@()joint_opening_profile_hydraulics(2e-10, ...
                [.001 .002],[2 3],struct('sampleWeights',[1 0])), ...
                'JointOpeningProfile:Weights');
        end

        function defaultMappingAndEmptyWarmTemperaturePreserveBase(testCase)
            cref=[2e-10;3e-10]; gap=[.001 .002;.002 .003];
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(cref,gap,[],struct());
            testCase.verifyEqual(C,cref,'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(2),'AbsTol',0);
            testCase.verifyEqual(Dt,zeros(2),'AbsTol',0);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,zeros(2,1),'AbsTol',0);
            testCase.verifyFalse(D.thermalProfileProvided);
        end

        function tinyPositiveConductanceIsNotFlooredOrLostToSquaredUnderflow(testCase)
            gap=[1e-85 2e-85 3e-85]; options=struct('openingTransmissionFraction',1, ...
                'baseMobilityFactor',0);
            [C,Dg,~,D]=joint_opening_profile_hydraulics(2e-10,gap,[2 3 4],options);
            testCase.verifyGreaterThan(C,0);
            testCase.verifyLessThan(C,1e-240);
            testCase.verifyGreaterThan(Dg,zeros(1,3));
            testCase.verifyTrue(all(isfinite(Dg)));
            testCase.verifyEqual(D.openConductance,C,'AbsTol',0);
        end

        function wholePathMaskBlocksFlowButPreservesPartiallyOpenInventory(testCase)
            options=struct('openingTransmissionFraction',.001, ...
                'baseMobilityFactor',.4,'baseMobilityDerivative_dT',[.1 -.2 .3]);
            gap=[.001 .002 .003]; T=[2 3 4];
            [~,~,~,unmasked]=joint_opening_profile_hydraulics(2e-10,gap,T,options);
            options.openPathMask=false;
            [C,Dg,Dt,D]=joint_opening_profile_hydraulics(2e-10,gap,T,options);
            testCase.verifyEqual(C,.4*2e-10,'AbsTol',0);
            testCase.verifyEqual(D.openConductance,0,'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(Dt,2e-10*[.1 -.2 .3],'AbsTol',1e-25);
            testCase.verifyEqual(D.sampleVolume_m2,unmasked.sampleVolume_m2,'AbsTol',0);
            testCase.verifyEqual(D.aperture_m,unmasked.aperture_m,'AbsTol',0);
            testCase.verifyGreaterThan(D.newslitGeometricVolume_m2,0);
        end

        function independentWholePathMasksPreserveJointInventoryAndIdentity(testCase)
            options=struct('openingTransmissionFraction',.001,'openPathMask',[false;true]);
            gap=[.001 .002;.002 .003]; T=4*ones(2);
            [C,Dg,~,D]=joint_opening_profile_hydraulics([2e-10;3e-10],gap,T,options);
            testCase.verifyEqual(C(1),2e-10,'AbsTol',0);
            testCase.verifyGreaterThan(C(2),3e-10);
            testCase.verifyEqual(Dg(1,:),zeros(1,2),'AbsTol',0);
            testCase.verifyGreaterThan(Dg(2,:),zeros(1,2));
            testCase.verifyGreaterThan(D.newslitGeometricVolume_m2,zeros(2,1));
        end
    end

    methods (Static,Access=private)
        function [cref,gap,T,options]=phaseFixture()
            cref=[2e-10;3e-10]; gap=[.001 .0018 .003;.0012 .0022 .0035];
            fraction=[.001;.0015]; gamma=.015;
            Tf=-2*273.15*gamma./(917*3.35e5*fraction.*gap);
            T=Tf+[.01 .025 .7;.02 .01 .5];
            options=struct('openingTransmissionFraction',fraction,'pathLength_m',[.45;.5], ...
                'sampleWeights',[.2 .3 .5],'freezeEnabled',true,'viscosityEnabled',true, ...
                'wallEnergyDifference_J_m2',gamma,'transitionWidth_C',.1);
        end

        function mu=viscosity(T)
            tau=(T+273.15)/300;
            mu=1e-6*(280.68*tau.^(-1.9)+511.45*tau.^(-7.7) ...
                +61.131*tau.^(-19.6)+.45903*tau.^(-40));
        end
    end
end
