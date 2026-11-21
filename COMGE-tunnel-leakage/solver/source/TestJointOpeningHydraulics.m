classdef TestJointOpeningHydraulics < matlab.unittest.TestCase
    %TESTJOINTOPENINGHYDRAULICS Independent slit flow/phase/derivative checks.
    % These validate ideal constitutive algebra, not a particular gasket.

    methods (Test)
        function defaultMappingPreservesBaseAndCannotInventLeaks(testCase)
            cref=[0;2e-10]; gap=[1e-3;1e-3];
            [C,Dg,Dt,D]=joint_opening_hydraulics(cref,gap,[],struct());
            testCase.verifyEqual(C,cref,'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(2,1),'AbsTol',0);
            testCase.verifySize(Dt,[2 0]);
            testCase.verifyEqual(D.openConductance,zeros(2,1),'AbsTol',0);
            testCase.verifyEqual(D.status{1},'sealed');
        end

        function sealedDefectRemainsExactlyClosedAtPositiveOpening(testCase)
            options=struct('openingTransmissionFraction',1,'freezeEnabled',true, ...
                'transitionWidth_C',.2,'wallEnergyDifference_J_m2',.029);
            [C,Dg,Dt,D]=joint_opening_hydraulics(0,.002,[-.01 2 5],options);
            testCase.verifyEqual(C,0,'AbsTol',0);
            testCase.verifyEqual(Dg,0,'AbsTol',0);
            testCase.verifyEqual(Dt,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(D.aperture_m,0,'AbsTol',0);
            testCase.verifyTrue(D.sealedMask);
        end

        function closedGapPreservesOriginalDefect(testCase)
            options=struct('openingTransmissionFraction',1);
            [C,Dg,~,D]=joint_opening_hydraulics([2e-10;3e-10],[-.001;0],[],options);
            testCase.verifyEqual(C,[2e-10;3e-10],'AbsTol',0);
            testCase.verifyEqual(Dg,zeros(2,1),'AbsTol',0);
            testCase.verifyEqual(D.openConductance,zeros(2,1),'AbsTol',0);
        end

        function warmSlitMatchesParallelPlateAnalyticalSolution(testCase)
            options=struct('openingTransmissionFraction',.002,'pathLength_m',.45);
            cref=2e-10; gap=.0008;
            mu=TestJointOpeningHydraulics.viscosity(20);
            aperture=gap*.002;
            expected=1000*9.81*aperture^3/(12*mu*.45);
            [C,Dg,~,D]=joint_opening_hydraulics(cref,gap,[],options);
            testCase.verifyEqual(C,cref+expected,'RelTol',1e-13);
            testCase.verifyEqual(D.openConductance,expected,'RelTol',1e-13);
            testCase.verifyEqual(Dg,3*expected/gap,'RelTol',1e-13);
            testCase.verifyEqual(D.openingContribution,aperture,'AbsTol',0);
            testCase.verifyGreaterThan(D.baseAperture,0);
        end

        function doublingGapMultipliesOnlyNewPassageByEight(testCase)
            options=struct('openingTransmissionFraction',.01);
            [~,~,~,D1]=joint_opening_hydraulics(1e-10,.001,[],options);
            [~,~,~,D2]=joint_opening_hydraulics(1e-10,.002,[],options);
            testCase.verifyEqual(D2.openConductance,8*D1.openConductance,'RelTol',1e-14);
            testCase.verifyEqual(D2.baseAperture,D1.baseAperture,'AbsTol',0);
        end

        function existingCapillaryMobilityIsPassedThroughSeparately(testCase)
            options=struct('openingTransmissionFraction',.001, ...
                'baseMobilityFactor',.3,'baseMobilityDerivative_dT',[.1 -.2 .3]);
            [C,Dg,Dt,D]=joint_opening_hydraulics(2e-10,.001,[1 2 3],options);
            testCase.verifyEqual(C,.3*2e-10+D.openConductance,'AbsTol',1e-24);
            testCase.verifyEqual(Dt,2e-10*[.1 -.2 .3],'AbsTol',1e-24);
            testCase.verifyGreaterThan(Dg,0);
        end

        function singleFrozenPathPointBlocksOpenSlitInSeries(testCase)
            options=struct('openingTransmissionFraction',.01,'freezeEnabled',true, ...
                'baseMobilityFactor',.4);
            [C,Dg,Dt,D]=joint_opening_hydraulics(2e-10,.001,[2 -.1 3],options);
            testCase.verifyEqual(C,.4*2e-10,'AbsTol',0);
            testCase.verifyEqual(D.openConductance,0,'AbsTol',0);
            testCase.verifyEqual(Dg,0,'AbsTol',0);
            testCase.verifyEqual(Dt,zeros(1,3),'AbsTol',0);
            testCase.verifyEqual(D.status{1},'open_slit_ice_blocked');
        end

        function slitGTUsesWallEnergyDifferenceAndFullAperture(testCase)
            gap=.001; fraction=.001; gamma=.015;
            options=struct('openingTransmissionFraction',fraction,'freezeEnabled',true, ...
                'wallEnergyDifference_J_m2',gamma);
            expectedTf=-2*273.15*gamma/(917*3.35e5*gap*fraction);
            [~,~,~,D]=joint_opening_hydraulics(2e-10,gap,[expectedTf+.001 1],options);
            testCase.verifyEqual(D.freezeTemperature_C,expectedTf,'RelTol',1e-14);
            testCase.verifyEqual(D.openMobilityFactor,1,'AbsTol',0);
            testCase.verifyEqual(D.conductingFractionBySample,ones(1,2),'AbsTol',0);
        end

        function temperatureDependentSeriesViscosityMatchesResistanceIntegral(testCase)
            T=[5 20 35]; weights=[.2 .3 .5];
            options=struct('openingTransmissionFraction',.001, ...
                'viscosityEnabled',true,'sampleWeights',weights);
            mu=TestJointOpeningHydraulics.viscosity(T);
            muRef=TestJointOpeningHydraulics.viscosity(20);
            expectedFactor=muRef/sum(weights.*mu);
            [~,~,~,D]=joint_opening_hydraulics(2e-10,.001,T,options);
            testCase.verifyEqual(D.openMobilityFactor,expectedFactor,'RelTol',1e-13);
        end

        function geometricalOpeningDerivativeMatchesFiniteDifference(testCase)
            options=struct('openingTransmissionFraction',[.001;.003]);
            cref=[2e-10;4e-10]; gap=[.001;.002]; v=[.2;-.3]; e=1e-8;
            [~,Dg]=joint_opening_hydraulics(cref,gap,[],options);
            p=joint_opening_hydraulics(cref,gap+e*v,[],options);
            m=joint_opening_hydraulics(cref,gap-e*v,[],options);
            testCase.verifyEqual(Dg.*v,(p-m)/(2*e),'RelTol',3e-8,'AbsTol',1e-15);
        end

        function regularizedSlitPhaseAndOpeningTangentsMatchFiniteDifference(testCase)
            options=struct('openingTransmissionFraction',.001,'freezeEnabled',true, ...
                'viscosityEnabled',true,'transitionWidth_C',.1, ...
                'wallEnergyDifference_J_m2',.015);
            cref=2e-10; gap=.001;
            Tf=-2*273.15*.015/(917*3.35e5*gap*.001);
            T=[Tf+.01 Tf+.025 .7]; vT=[.13 -.07 .04]; eT=1e-6; eGap=1e-8;
            [~,Dg,Dt,D]=joint_opening_hydraulics(cref,gap,T,options);
            pT=joint_opening_hydraulics(cref,gap,T+eT*vT,options);
            mT=joint_opening_hydraulics(cref,gap,T-eT*vT,options);
            pg=joint_opening_hydraulics(cref,gap+eGap,T,options);
            mg=joint_opening_hydraulics(cref,gap-eGap,T,options);
            testCase.verifyEqual(Dt*vT',(pT-mT)/(2*eT),'RelTol',2e-7,'AbsTol',1e-18);
            testCase.verifyEqual(Dg,(pg-mg)/(2*eGap),'RelTol',3e-7,'AbsTol',1e-15);
            testCase.verifyGreaterThan(D.openMobilityFactor,0);
            testCase.verifyLessThan(D.openMobilityFactor,1);
        end

        function missingThermalPathAndInvalidMappingAreRejected(testCase)
            testCase.verifyError(@()joint_opening_hydraulics(1e-10,.001,[], ...
                struct('freezeEnabled',true)),'JointOpeningHydraulics:MissingTemperature');
            testCase.verifyError(@()joint_opening_hydraulics(1e-10,.001,[], ...
                struct('openingTransmissionFraction',1.1)),'JointOpeningHydraulics:JointParameter');
        end

        function channelVolumesPartitionAndDifferentiateExactly(testCase)
            options=struct('openingTransmissionFraction',.001,'freezeEnabled',true, ...
                'transitionWidth_C',.1,'wallEnergyDifference_J_m2',.015, ...
                'sampleWeights',[.2 .3 .5]);
            gap=.001; Tf=-2*273.15*.015/(917*3.35e5*gap*.001);
            T=[Tf+.01 Tf+.025 .7]; vT=[.13 -.07 .04]; eT=1e-6; eg=1e-8;
            [~,~,~,D]=joint_opening_hydraulics(2e-10,gap,T,options);
            [~,~,~,pT]=joint_opening_hydraulics(2e-10,gap,T+eT*vT,options);
            [~,~,~,mT]=joint_opening_hydraulics(2e-10,gap,T-eT*vT,options);
            [~,~,~,pg]=joint_opening_hydraulics(2e-10,gap+eg,T,options);
            [~,~,~,mg]=joint_opening_hydraulics(2e-10,gap-eg,T,options);
            testCase.verifyEqual(D.newslitGeometricVolume_m2,.45*.001*gap,'RelTol',1e-14);
            testCase.verifyEqual(D.liquidVolumeBySample_m2+D.iceVolumeBySample_m2, ...
                D.sampleVolume_m2,'AbsTol',1e-21);
            testCase.verifyEqual(sum(D.sampleVolume_m2),D.newslitGeometricVolume_m2, ...
                'AbsTol',1e-21);
            testCase.verifyEqual(D.dLiquidVolume_dGap_m+D.dIceVolume_dGap_m, ...
                D.dVolume_dGap_m*D.sampleWeights,'AbsTol',1e-18);
            testCase.verifyEqual(D.dLiquidVolume_dGap_m, ...
                (pg.liquidVolumeBySample_m2-mg.liquidVolumeBySample_m2)/(2*eg), ...
                'RelTol',3e-7,'AbsTol',1e-12);
            testCase.verifyEqual(D.dLiquidVolume_dTemperature_m2_K.*vT, ...
                (pT.liquidVolumeBySample_m2-mT.liquidVolumeBySample_m2)/(2*eT), ...
                'RelTol',2e-7,'AbsTol',1e-15);
        end
    end

    methods (Static,Access=private)
        function mu=viscosity(T)
            tau=(T+273.15)/300;
            mu=1e-6*(280.68*tau.^(-1.9)+511.45*tau.^(-7.7) ...
                +61.131*tau.^(-19.6)+.45903*tau.^(-40));
        end
    end
end
