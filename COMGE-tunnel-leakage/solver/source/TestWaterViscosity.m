classdef TestWaterViscosity < matlab.unittest.TestCase
    methods (Test)
        function reproducesAllElevenIapwsTable4Points(tc)
            T=[298.15,298.15,373.15,433.15,433.15,873.15,873.15,873.15,1173.15,1173.15,1173.15]';
            rho=[998,1200,1000,1,1000,1,100,600,1,100,400]';
            expected=1e-6*[889.735100,1437.649467,307.883622,14.538324,217.685358,32.619287,35.802262,77.430195,44.217245,47.640433,64.154608]';
            tc.verifyEqual(water_viscosity_iapws2008(T,rho),expected,'AbsTol',5.1e-13);
        end
        function analyticDerivativesMatchIndependentFiniteDifferences(tc)
            T=[269.15;279.15;293.15]; rho=[990;1000;1020]; h=.001;
            [~,a,b]=water_viscosity_iapws2008(T,rho);
            fa=(water_viscosity_iapws2008(T+h,rho)-water_viscosity_iapws2008(T-h,rho))/(2*h);
            fb=(water_viscosity_iapws2008(T,rho+h)-water_viscosity_iapws2008(T,rho-h))/(2*h);
            tc.verifyLessThan(norm(a-fa)/norm(a),1e-7); tc.verifyLessThan(norm(b-fb)/norm(b),1e-7);
        end
        function nearCriticalAndUnsupportedTemperatureAreRejected(tc)
            tc.verifyError(@()water_viscosity_iapws2008(647,322),'WaterViscosity:Domain');
            tc.verifyError(@()water_viscosity_iapws2008(240,1000),'WaterViscosity:Domain');
        end
    end
end
