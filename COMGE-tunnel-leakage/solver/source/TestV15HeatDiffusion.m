classdef TestV15HeatDiffusion < matlab.unittest.TestCase
    % Actual solver verification against a continuous analytic heat solution.
    properties (Access=private)
        study
    end
    methods (TestClassSetup)
        function solveOneSharedBenchmark(testCase)
            % One study serves every test; no repeated transient solves.
            testCase.study=run_v15_heat_diffusion_benchmark(struct('verbose',false));
        end
    end
    methods (Test)
        function warmDiffusionConvergesUnderMeshRefinement(testCase)
            r=testCase.study.comparison;
            testCase.verifyLessThan(diff(r.L2_C),zeros(height(r)-1,1));
            testCase.verifyGreaterThan(r.observedOrder(end),1.6);
            testCase.verifyLessThan(r.L2_C(end),3e-3);
        end
        function temporalContaminationIsSmall(testCase)
            t=testCase.study.timeRefinement;
            testCase.verifyTrue(t.performed);
            testCase.verifyLessThan(t.estimatedTemporalFractionOfFineL2,.10);
            testCase.verifyEqual(t.fine.metrics.dtSeconds, ...
                testCase.study.comparison.dtSeconds(end)/2);
        end
        function hydrostaticPatchDoesNotCreateAdvection(testCase)
            r=testCase.study.comparison;
            testCase.verifyLessThan(r.maxDarcyVelocity_m_s,1e-13*ones(height(r),1));
            testCase.verifyLessThan(r.maxSubfaceFlow_m2_s,1e-13*ones(height(r),1));
            testCase.verifyLessThan(r.maxHydrostaticHeadError_m,1e-7*ones(height(r),1));
            for j=1:numel(testCase.study.caseResults)
                testCase.verifyTrue(testCase.study.caseResults{j}.out.parameters.advection);
            end
        end
        function boundaryAndPhaseRemainTheAnalyticWarmCase(testCase)
            r=testCase.study.comparison;
            testCase.verifyLessThan(r.maxTopTemperatureError_C,1e-10*ones(height(r),1));
            testCase.verifyEqual(r.maxIceVolumeFraction,zeros(height(r),1),'AbsTol',0);
            testCase.verifyLessThan(r.relativeElasticHeatCapacity,1e-14*ones(height(r),1));
            testCase.verifyEqual(r.rejectedSteps,zeros(height(r),1));
        end
        function waterAndEnergyLedgersClose(testCase)
            r=testCase.study.comparison;
            testCase.verifyLessThan(abs(r.cumulativeWaterDefect_m2),1e-10*ones(height(r),1));
            testCase.verifyLessThan(r.relativeEnergyBudgetDefect,1e-8*ones(height(r),1));
        end
    end
end
