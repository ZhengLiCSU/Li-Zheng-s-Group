classdef TestGasketLiteratureReference < matlab.unittest.TestCase
    % Compare published tables, not a refitted synthetic reference curve.
    methods (Test)
        function agingFormulaMatchesPublishedRoundedTable5(testCase)
            result = zhang2020_gasket_reference([0 5 10 30 50 75 100],7);
            expected = [1 0.748 0.679 0.542 0.469 0.408 0.364];
            testCase.verifyEqual(result.materialRetention,expected,'AbsTol',0.00051);
        end
        function waterproofFormulaMatchesPublishedTable9(testCase)
            result = zhang2020_gasket_reference(100,[7 9 10]);
            testCase.verifyEqual(result.waterproofPressureRetention,[0.237 0.312 0.358], ...
                'AbsTol',0.00051);
        end
        function unagedRatioIsUnityAtAllTestedOpenings(testCase)
            result = zhang2020_gasket_reference(0,[7 9 10]);
            testCase.verifyEqual(result.waterproofPressureRetention,ones(1,3),'AbsTol',0);
        end
        function unsupportedOpeningIsRejected(testCase)
            testCase.verifyError(@()zhang2020_gasket_reference(20,0),'GasketReference:Opening');
        end
        function unsupportedAgeIsRejected(testCase)
            testCase.verifyError(@()zhang2020_gasket_reference(101,7),'GasketReference:Age');
        end
        function temperatureAndScopeAreExplicit(testCase)
            result = zhang2020_gasket_reference(20,7);
            testCase.verifyEqual(result.serviceTemperature_K,293);
            testCase.verifyTrue(contains(result.scope,'no freeze-thaw calibration'));
        end
    end
end
