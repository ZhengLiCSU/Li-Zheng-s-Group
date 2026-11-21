classdef verify_joint_freeze_law < matlab.unittest.TestCase
    %VERIFY_JOINT_FREEZE_LAW Independent limits and tangent checks.
    % Run with results=runtests('verify_joint_freeze_law.m');
    %          assertSuccess(results);
    % No result here is a calibrated tunnel prediction or an expected
    % winter shut-off date. These checks exercise public physical limits.

    properties (TestParameter)
        derivativeCase = struct( ...
            'uniformWarm',20*ones(2,8), ...
            'uniformTransition',-.05*ones(2,8), ...
            'mixedPath',[-.10 -.06 -.04 -.02 0 1 5 20; ...
                         10 5 2 .1 -.015 -.05 -.1 -.3]);
    end

    methods (Test)
        function disabledIsExactIdentity(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=false;
            profile=[-5 -1 -.05 0 1 2 10 20; ...
                     20 10 2 1 0 -.05 -1 -5];
            [factor,derivative]=joint_freeze_conductance(profile,settings);
            testCase.verifyEqual(factor,ones(2,1),'AbsTol',0);
            testCase.verifyEqual(derivative,zeros(size(profile)),'AbsTol',0);
        end

        function referenceWarmTemperatureHasUnitFactor(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            profile=settings.referenceTemperature_C*ones(2,8);
            [factor,derivative,diagnostic]=joint_freeze_conductance(profile,settings);
            testCase.verifyEqual(factor,ones(2,1),'AbsTol',1e-13);
            testCase.verifyTrue(all(isfinite(derivative(:))));
            testCase.verifyEqual(diagnostic.iceConductanceFraction, ...
                ones(size(profile)),'AbsTol',1e-13);
        end

        function veryColdAndTransitionStatesStayFinite(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            temperature=[-20;-5;-1;-.2;-.051;-.05;-.049;-.001;0;.001;20];
            profile=repmat(temperature,1,8);
            [factor,derivative,diagnostic]=joint_freeze_conductance(profile,settings);
            testCase.verifyTrue(all(isfinite(factor(:))) && all(factor(:)>=0));
            testCase.verifyTrue(all(isfinite(derivative(:))));
            testCase.verifyTrue(all(isfinite(diagnostic.mobility(:))) && ...
                all(diagnostic.mobility(:)>=0));
            testCase.verifyTrue(all(diagnostic.iceConductanceFraction(:)>=0) && ...
                all(diagnostic.iceConductanceFraction(:)<=1));
        end

        function coldBottleneckLimitsPositiveMeanPath(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            settings.viscosityEnabled=false;
            profile=[-.1 2 2 2 2 2 2 2];
            [factor,~,diagnostic]=joint_freeze_conductance(profile,settings);
            meanFactor=joint_freeze_conductance(mean(profile)*ones(1,8),settings);
            testCase.verifyGreaterThan(mean(profile),0);
            testCase.verifyLessThan(factor,.95);
            testCase.verifyLessThan(factor,meanFactor);
            testCase.verifyLessThan(factor,mean(diagnostic.mobility,2));
            testCase.verifyEqual(factor,1./mean(1./diagnostic.mobility,2), ...
                'AbsTol',1e-12);
        end

        function uniformWarmingIsMonotone(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            temperature=linspace(-2,20,501)';
            [factor,derivative]=joint_freeze_conductance(repmat(temperature,1,8),settings);
            testCase.verifyGreaterThanOrEqual(min(diff(factor)),-1e-12);
            testCase.verifyGreaterThanOrEqual(min(derivative(:)),-1e-12);
            testCase.verifyGreaterThan(factor(end),factor(1));
        end

        function reversibleLawReturnsOriginalWarmConductance(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            warm=5*ones(1,8);
            before=joint_freeze_conductance(warm,settings);
            cold=joint_freeze_conductance(-.1*ones(1,8),settings);
            after=joint_freeze_conductance(warm,settings);
            testCase.verifyLessThan(cold,before);
            testCase.verifyEqual(after,before,'AbsTol',0);
        end

        function everyPathTemperatureDerivativeMatchesFD(testCase,derivativeCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            [factor,derivative]=joint_freeze_conductance(derivativeCase,settings);
            [fd,analytic]=verify_joint_freeze_law.finiteDifference( ...
                derivativeCase,settings,derivative);
            relativeError=norm(fd-analytic,inf)/max(norm(fd,inf),1e-10);
            testCase.verifySize(factor,[size(derivativeCase,1),1]);
            testCase.verifySize(derivative,size(derivativeCase));
            testCase.verifyLessThan(relativeError,2e-5);
        end

        function perJointRadiiEqualIndependentScalarCalls(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            settings.viscosityEnabled=false;
            settings.flowWeightedRadius_m=[.5e-6;2e-6];
            settings.logRadiusStd=[.5;.9];
            profile=-.05*ones(2,8);
            [factor,derivative]=joint_freeze_conductance(profile,settings);
            first=settings;
            first.flowWeightedRadius_m=.5e-6; first.logRadiusStd=.5;
            second=settings;
            second.flowWeightedRadius_m=2e-6; second.logRadiusStd=.9;
            [f1,d1]=joint_freeze_conductance(profile(1,:),first);
            [f2,d2]=joint_freeze_conductance(profile(2,:),second);
            testCase.verifyEqual(factor,[f1;f2],'AbsTol',1e-13);
            testCase.verifyEqual(derivative,[d1;d2],'AbsTol',1e-12);
        end

        function viscosityAndIceCanBeSeparated(testCase)
            settings=joint_freeze_defaults();
            settings.enabled=true;
            profile=-.05*ones(1,8);
            full=joint_freeze_conductance(profile,settings);
            viscosity=settings; viscosity.iceEnabled=false;
            ice=settings; ice.viscosityEnabled=false;
            neutral=ice; neutral.iceEnabled=false;
            [viscosityOnly,~,viscosityDiagnostic]=joint_freeze_conductance(profile,viscosity);
            iceOnly=joint_freeze_conductance(profile,ice);
            [neutralFactor,neutralDerivative]=joint_freeze_conductance(profile,neutral);
            testCase.verifyGreaterThan(viscosityOnly,0);
            testCase.verifyLessThan(viscosityOnly,1);
            testCase.verifyGreaterThan(iceOnly,0);
            testCase.verifyLessThan(iceOnly,1);
            testCase.verifyEqual(viscosityDiagnostic.iceConductanceFraction, ...
                ones(size(profile)),'AbsTol',1e-13);
            testCase.verifyEqual(full,viscosityOnly*iceOnly,'AbsTol',1e-12);
            testCase.verifyEqual(neutralFactor,1,'AbsTol',0);
            testCase.verifyEqual(neutralDerivative,zeros(size(profile)),'AbsTol',0);
        end
    end

    methods (Static,Access=private)
        function [fd,analytic]=finiteDifference(profile,settings,derivative)
            % Uncompressed differences detect both local derivative errors
            % and accidental cross-joint coupling, independently of G.jnz.
            rows=size(profile,1); columns=size(profile,2);
            fd=zeros(rows,numel(profile));
            for index=1:numel(profile)
                increment=1e-6*max(1,abs(profile(index)));
                plus=profile; minus=profile;
                plus(index)=plus(index)+increment;
                minus(index)=minus(index)-increment;
                fp=joint_freeze_conductance(plus,settings);
                fm=joint_freeze_conductance(minus,settings);
                fd(:,index)=(fp-fm)/(2*increment);
            end
            analytic=full(sparse(repmat((1:rows)',columns,1), ...
                (1:numel(profile))',derivative(:),rows,numel(profile)));
        end
    end
end
