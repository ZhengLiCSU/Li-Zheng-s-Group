classdef TestSoilSCVQuadrature < matlab.unittest.TestCase
    properties (TestParameter)
        order = {2,5}
    end
    methods (Test)
        function constantVolumes(testCase,order)
            Q=build_soil_scv_quadrature([0 0;2 0;.3 1],[1 2 3],order);
            testCase.verifyEqual(Q.volume,ones(3,1)/3,'AbsTol',2e-14);
            testCase.verifyEqual(full(sum(Q.P,2)),ones(size(Q.P,1),1),'AbsTol',1e-14);
            testCase.verifyEqual(full(Q.elementAverage*ones(size(Q.P,1),1)),1,'AbsTol',1e-14);
            testCase.verifyGreaterThan(Q.weight,0);
        end
        function affineInterpolation(testCase,order)
            xy=[0 0;2 0;.3 1]; Q=build_soil_scv_quadrature(xy,[1 2 3],order);
            testCase.verifyEqual(Q.P*(3+2*xy(:,1)-xy(:,2)), ...
                3+2*Q.coordinates(:,1)-Q.coordinates(:,2),'AbsTol',2e-14);
        end
        function quadraticIntegral(testCase,order)
            Q=build_soil_scv_quadrature([0 0;1 0;0 1],[1 2 3],order);
            testCase.verifyEqual(sum(Q.weight.*Q.coordinates(:,1).^2),1/12,'AbsTol',1e-14);
            testCase.verifyEqual(sum(Q.weight.*prod(Q.coordinates,2)),1/24,'AbsTol',1e-14);
        end
        function physicalVolumeIsNotGalerkinWeight(testCase)
            Q=build_soil_scv_quadrature([0 0;1 0;0 1],[1 2 3],5);
            integral=Q.deposit*(Q.P*[0;1;0]);
            testCase.verifyEqual(sum(integral),1/6,'AbsTol',1e-14);
            testCase.verifyNotEqual(integral,[1;2;1]/24);
            testCase.verifyGreaterThan(integral(2),integral(1));
        end
        function reversedOrientation(testCase,order)
            xy=[0 0;1 0;0 1]; Q1=build_soil_scv_quadrature(xy,[1 2 3],order);
            Q2=build_soil_scv_quadrature(xy,[1 3 2],order);
            u=[.2;2;-.3];
            testCase.verifyEqual(Q1.deposit*(Q1.P*u).^2,Q2.deposit*(Q2.P*u).^2,'AbsTol',1e-14);
        end
        function invalidOrder(testCase)
            testCase.verifyError(@()build_soil_scv_quadrature([0 0;1 0;0 1],[1 2 3],3),'SCV:order');
        end
        function subdivisionPreservesOwnerIntegral(testCase)
            xy=[0 0;1 0;0 1]; u=[.2;2;-.3];
            Q0=build_soil_scv_quadrature(xy,[1 2 3],5,0);
            Q1=build_soil_scv_quadrature(xy,[1 2 3],5,1);
            Q2=build_soil_scv_quadrature(xy,[1 2 3],5,2);
            testCase.verifyEqual(Q1.deposit*(Q1.P*u).^2,Q0.deposit*(Q0.P*u).^2,'AbsTol',1e-14);
            testCase.verifyEqual(Q2.deposit*(Q2.P*u).^2,Q0.deposit*(Q0.P*u).^2,'AbsTol',1e-14);
            testCase.verifyEqual(size(Q2.P,1),16*size(Q0.P,1));
        end
    end
end
