classdef TestJointChannelInventory < matlab.unittest.TestCase
    %TESTJOINTCHANNELINVENTORY Independent capillary storage and heat checks.
    % Conductance/area is checked from a known number of lognormal tubes,
    % and energy/head tangents by finite differences of public interfaces.
    properties (TestParameter)
        invalidRadius=struct('zero',0,'negative',-1e-6);
        invalidStd=struct('zero',0,'negative',-.7);
        invalidTortuosity=struct('zero',0,'negative',-1,'belowUnity',.9);
    end

    methods (Test)
        function flowWeightedToVolumeWeightedFreezing(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            beta=2*freeze.interfacialEnergy_J_m2*freeze.meltingTemperature_K/ ...
                (freeze.iceDensity_kg_m3*freeze.latentHeat_J_kg);
            T=freeze.bulkFreezingTemperature_C-beta/freeze.flowWeightedRadius_m;
            [~,~,info]=TestJointChannelInventory.evaluate(T,freeze,storage,k);
            factor=joint_freeze_conductance(T,freeze);
            expectedVolume=.5*erfc(-2*freeze.logRadiusStd/sqrt(2));
            testCase.verifyEqual(factor,.5,'AbsTol',2e-14);
            testCase.verifyEqual(info.liquidVolumeFraction,expectedVolume,'AbsTol',2e-14);
            testCase.verifyGreaterThan(info.liquidVolumeFraction,factor);
            testCase.verifyEqual(info.volumeWeightedRadius_m, ...
                freeze.flowWeightedRadius_m*exp(-2*freeze.logRadiusStd^2), ...
                'AbsTol',1e-18);
        end

        function knownLognormalTubePopulationMatchesAreaAndVolume(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            numberPerMetre=1e7; medianRadius=1e-6; sigma=.7;
            freeze.flowWeightedRadius_m=medianRadius*exp(4*sigma^2);
            freeze.logRadiusStd=sigma; storage.tortuosity=1.4;
            pathLength=k.L*storage.tortuosity;
            mu=water_viscosity_iapws(freeze.referenceTemperature_C);
            moment2=medianRadius^2*exp(2*sigma^2);
            moment4=medianRadius^4*exp(8*sigma^2);
            k.c=k.rho*k.g*numberPerMetre*pi*moment4/(8*mu*pathLength);
            [~,~,info]=TestJointChannelInventory.evaluate([20 20],freeze,storage,k);
            expectedArea=numberPerMetre*pi*moment2;
            testCase.verifyEqual(info.connectedArea_m2_per_m,expectedArea, ...
                'AbsTol',1e-16,'RelTol',1e-13);
            testCase.verifyEqual(info.waterVolume_m2,expectedArea*pathLength, ...
                'AbsTol',1e-16,'RelTol',1e-13);
            testCase.verifyEqual(info.pathLength_m,pathLength,'AbsTol',1e-15);
        end

        function fixedConductanceVolumeScalesWithSquaredPathRatio(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            [~,~,a]=TestJointChannelInventory.evaluate(5,freeze,storage,k);
            storage.tortuosity=2*storage.tortuosity;
            [~,~,b]=TestJointChannelInventory.evaluate(5,freeze,storage,k);
            testCase.verifyEqual(b.connectedArea_m2_per_m,2*a.connectedArea_m2_per_m, ...
                'AbsTol',1e-16,'RelTol',1e-13);
            testCase.verifyEqual(b.waterVolume_m2,4*a.waterVolume_m2, ...
                'AbsTol',1e-16,'RelTol',1e-13);
        end

        function warmEnergyIncludesLatentOffsetAndDisplacedConcrete(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            T=[5 12 20];
            [E,JT,info]=TestJointChannelInventory.evaluate(T,freeze,storage,k);
            sampleVolume=info.waterVolume_m2/3;
            expectedE=sampleVolume*(k.rho*k.latent+(k.rho*k.cw-k.rhoCcc)*T);
            expectedJ=sampleVolume*(k.rho*k.cw-k.rhoCcc)*ones(1,3);
            testCase.verifyEqual(E,expectedE,'AbsTol',1e-9,'RelTol',1e-13);
            testCase.verifyEqual(JT,expectedJ,'AbsTol',1e-10,'RelTol',1e-13);
            testCase.verifyEqual(info.liquidVolumeFraction,ones(1,3),'AbsTol',0);
        end

        function energyDerivativeMatchesFiniteDifferenceThroughFreezing(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            freeze.flowWeightedRadius_m=[1e-6;2e-6];
            freeze.logRadiusStd=[.7;1.1]; storage.tortuosity=[1.2;1.8];
            T=[-.02 -.08 -.3 4;-.04 -.1 -1 12];
            direction=[.2 -.3 .4 -.1;-.3 .1 -.2 .3]; epsilon=1e-7;
            [~,JT]=TestJointChannelInventory.evaluate(T,freeze,storage,k);
            Ep=TestJointChannelInventory.evaluate(T+epsilon*direction,freeze,storage,k);
            Em=TestJointChannelInventory.evaluate(T-epsilon*direction,freeze,storage,k);
            testCase.verifyEqual(JT.*direction,(Ep-Em)/(2*epsilon), ...
                'AbsTol',2e-5,'RelTol',2e-6);
        end

        function uniformProfileInventoryDoesNotDependOnSampleCount(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            [E1,J1]=TestJointChannelInventory.evaluate(-.06,freeze,storage,k);
            [E8,J8]=TestJointChannelInventory.evaluate(-.06*ones(1,8),freeze,storage,k);
            testCase.verifyEqual(sum(E8),E1,'AbsTol',1e-10,'RelTol',1e-13);
            testCase.verifyEqual(sum(J8),J1,'AbsTol',1e-10,'RelTol',1e-13);
        end

        function explicitGravityOverridesDefaultAndPreservesScaling(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture();
            [Edefault,Jdefault,idefault]=joint_channel_inventory([-.04 5],freeze,storage, ...
                k.c,k.L,k.rho,k.latent,k.cw,k.ci,k.rhoCcc);
            k.g=2*9.81;
            [Echanged,Jchanged,ichanged]=TestJointChannelInventory.evaluate( ...
                [-.04 5],freeze,storage,k);
            testCase.verifyEqual(Echanged,Edefault/2,'AbsTol',1e-10,'RelTol',1e-13);
            testCase.verifyEqual(Jchanged,Jdefault/2,'AbsTol',1e-10,'RelTol',1e-13);
            testCase.verifyEqual(ichanged.waterVolume_m2,idefault.waterVolume_m2/2, ...
                'AbsTol',1e-16,'RelTol',1e-13);
        end

        function disabledStorageReturnsExactZeroInventory(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); storage.enabled=false;
            [E,JT,info]=TestJointChannelInventory.evaluate([-.1 5],freeze,storage,k);
            testCase.verifyEqual(E,zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(JT,zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(info.waterVolume_m2,0,'AbsTol',0);
            testCase.verifyEqual(info.connectedArea_m2_per_m,0,'AbsTol',0);
        end

        function zeroConductanceReturnsExactZeroWithoutZeroOverZero(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); k.c=[0;0];
            [E,JT,info]=TestJointChannelInventory.evaluate( ...
                [-20 -.05 0 20;-10 -.1 3 10],freeze,storage,k);
            testCase.verifyEqual(E,zeros(2,4),'AbsTol',0);
            testCase.verifyEqual(JT,zeros(2,4),'AbsTol',0);
            testCase.verifyEqual(info.waterVolume_m2,zeros(2,1),'AbsTol',0);
            testCase.verifyTrue(all(isfinite(info.temperatureDerivative_per_C(:))));
        end

        function disabledFreezingRetainsLiquidInventory(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); freeze.enabled=false;
            [~,JT,info]=TestJointChannelInventory.evaluate([-10 -.1],freeze,storage,k);
            testCase.verifyEqual(info.liquidVolumeFraction,ones(1,2),'AbsTol',0);
            testCase.verifyEqual(info.temperatureDerivative_per_C,zeros(1,2),'AbsTol',0);
            testCase.verifyEqual(JT,info.waterVolume_m2/2* ...
                (k.rho*k.cw-k.rhoCcc)*ones(1,2),'AbsTol',1e-11,'RelTol',1e-13);
        end

        function zeroConductanceDoesNotHideInvalidRadius(testCase,invalidRadius)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); k.c=0;
            freeze.flowWeightedRadius_m=invalidRadius;
            testCase.verifyError(@()TestJointChannelInventory.evaluate( ...
                -1,freeze,storage,k),'joint_channel_inventory:InvalidRadius');
        end

        function disabledStorageDoesNotHideInvalidRadiusStd(testCase,invalidStd)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); storage.enabled=false;
            freeze.logRadiusStd=invalidStd;
            testCase.verifyError(@()TestJointChannelInventory.evaluate( ...
                -1,freeze,storage,k),'joint_channel_inventory:InvalidLogRadiusStd');
        end

        function zeroConductanceDoesNotHideInvalidPathRatio(testCase,invalidTortuosity)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); k.c=0;
            storage.tortuosity=invalidTortuosity;
            testCase.verifyError(@()TestJointChannelInventory.evaluate( ...
                -1,freeze,storage,k),'joint_channel_inventory:InvalidTortuosity');
        end

        function negativeConductanceIsRejected(testCase)
            [freeze,storage,k]=TestJointChannelInventory.fixture(); k.c=-1e-10;
            testCase.verifyError(@()TestJointChannelInventory.evaluate( ...
                5,freeze,storage,k),'joint_channel_inventory:InvalidConductance');
        end

        function consistentInletHeatHeadTangentMatchesFiniteDifference(testCase)
            [T,S,paths,M,h]=TestJointChannelInventory.heatFixture();
            flow=sparse((M*h)'); slope=sparse(sum(M,1));
            inletDerivative=sparse((M'*(S*T))'); direction=[.2;-.4]; epsilon=1e-6;
            [~,~,J]=joint_channel_heat(T,S,paths,flow,flow,slope,sparse(1,4),1,inletDerivative);
            hp=joint_channel_heat(T,S,paths,sparse((M*(h+epsilon*direction))'),[],[],[],1);
            hm=joint_channel_heat(T,S,paths,sparse((M*(h-epsilon*direction))'),[],[],[],1);
            [~,~,legacyJ]=joint_channel_heat(T,S,paths,flow,flow,slope,sparse(1,4),1);
            testCase.verifyEqual(full(J*direction),full((hp-hm)/(2*epsilon)),'AbsTol',2e-9);
            testCase.verifyGreaterThan(norm(J-legacyJ,inf),.1);
        end

        function explicitLegacyInletDerivativePreservesEightInputAPI(testCase)
            [T,S,paths,M,h]=TestJointChannelInventory.heatFixture();
            flow=sparse((M*h)'); slope=sparse(sum(M,1));
            legacyDerivative=slope*spdiags(S*T,0,2,2);
            [oldHeat,oldT,oldH,oldExport]=joint_channel_heat( ...
                T,S,paths,flow,flow,slope,sparse(1,4),1);
            [newHeat,newT,newH,newExport]=joint_channel_heat( ...
                T,S,paths,flow,flow,slope,sparse(1,4),1,legacyDerivative);
            testCase.verifyEqual(newHeat,oldHeat,'AbsTol',0);
            testCase.verifyEqual(newT,oldT,'AbsTol',0);
            testCase.verifyEqual(newH,oldH,'AbsTol',0);
            testCase.verifyEqual(newExport,oldExport,'AbsTol',0);
        end

        function legacyHeatSelfTestStillPasses(testCase)
            report=joint_channel_heat('selftest');
            testCase.verifyTrue(report.passed);
            testCase.verifyLessThan(report.headDirectionalRelativeError,1e-7);
        end
    end

    methods (Static,Access=private)
        function [freeze,storage,k]=fixture()
            freeze=joint_freeze_defaults(struct('enabled',true, ...
                'viscosityEnabled',false,'flowWeightedRadius_m',1e-6,'logRadiusStd',.7));
            storage=struct('enabled',true,'tortuosity',1);
            k=struct('c',2e-10,'L',.45,'rho',1000,'latent',3.34e5, ...
                'cw',4180,'ci',2100,'rhoCcc',2500*900,'g',9.81);
        end

        function [E,JT,info]=evaluate(T,freeze,storage,k)
            [E,JT,info]=joint_channel_inventory(T,freeze,storage, ...
                k.c,k.L,k.rho,k.latent,k.cw,k.ci,k.rhoCcc,k.g);
        end

        function [T,S,paths,M,h]=heatFixture()
            T=[8;4;3;1]; S=sparse(1:2,1:2,1,2,4);
            paths=struct('nj',1,'nq',2,'interpolation',sparse([1 2],[3 4],1,2,4));
            M=[2 1;1 2]/6; h=[2;6];
        end
    end
end
