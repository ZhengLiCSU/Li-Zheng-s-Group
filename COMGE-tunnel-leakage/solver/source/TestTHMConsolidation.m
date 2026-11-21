classdef TestTHMConsolidation < matlab.unittest.TestCase
    % Analytic coupling verification; this is not a complete tunnel test.
    properties (Access=private)
        Study
    end
    methods (TestClassSetup)
        function runOneSharedBenchmark(testCase)
            testCase.Study=run_thm_consolidation_benchmark( ...
                struct('saveArtifacts',false,'verbose',false));
        end
    end
    methods (Test)
        function constrainedModulusAndDiffusivityMatchIndependentTheory(testCase)
            s=testCase.Study; E=80e6; nu=.35; biot=.65; rhoG=9810;
            modulus=E*(1-nu)/((1+nu)*(1-2*nu));
            storage=1e-5+rhoG*biot^2/modulus;
            testCase.verifyEqual(s.theory.constrainedModulus_Pa,modulus,'RelTol',1e-14);
            testCase.verifyEqual(s.theory.effectiveStorage_per_m,storage,'RelTol',1e-14);
            testCase.verifyEqual(s.theory.diffusivity_m2_s,1e-7/storage,'RelTol',1e-14);
            testCase.verifyGreaterThan(s.theory.poroelasticStorage_per_m,1e-5);
        end

        function spatialHeadErrorsDecreaseAtApproximatelySecondOrder(testCase)
            r=testCase.Study.comparison;
            testCase.verifyLessThan(diff(r.rmsHeadError_m),zeros(height(r)-1,1));
            testCase.verifyGreaterThan(r.observedOrder(end),1.6);
            testCase.verifyLessThan(r.observedOrder(end),2.6);
            testCase.verifyLessThan(r.relativeL2HeadError(end),.02);
        end

        function temporalContaminationIsMeasuredAndSmall(testCase)
            t=testCase.Study.timeRefinement;
            testCase.verifyTrue(t.performed);
            testCase.verifyGreaterThan(t.rmsHeadChange_m,0);
            testCase.verifyLessThan(t.estimatedTemporalFractionOfFineL2,.2);
            testCase.verifyEqual(t.fine.metrics.dtSeconds,5,'AbsTol',1e-12);
        end

        function waterStorageAndDarcyReactionsCloseLocallyAndGlobally(testCase)
            r=testCase.Study.comparison;
            testCase.verifyLessThan(r.relativeCumulativeWaterDefect,1e-8*ones(height(r),1));
            testCase.verifyLessThan(r.relativeSumAbsWaterDefect,1e-7*ones(height(r),1));
            testCase.verifyLessThan(r.maxLocalWaterDefect_m2,1e-12*ones(height(r),1));
            testCase.verifyEqual(r.totalStorageChange_m2,-r.totalDrainedWater_m2, ...
                'AbsTol',1e-12,'RelTol',1e-8);
            testCase.verifyGreaterThan(r.totalDrainedWater_m2,zeros(height(r),1));
        end

        function mechanicalLoadStorageAdjointAndInstalledReferenceAreConsistent(testCase)
            r=testCase.Study.comparison;
            testCase.verifyLessThan(r.storageLoadAdjointRelativeError,1e-12*ones(height(r),1));
            testCase.verifyLessThan(r.installedReferenceRelativeResidual,1e-11*ones(height(r),1));
            testCase.verifyLessThan(r.relativeMechanicalResidual,1e-10*ones(height(r),1));
            testCase.verifyLessThan(r.maxAbsIncrementalVolumeStrain,1e-4*ones(height(r),1));
        end

        function pressureRemainsSaturatedAndPressureLossCausesSettlement(testCase)
            r=testCase.Study.comparison;
            testCase.verifyGreaterThan(r.minimumAbsolutePressureHead_m,ones(height(r),1));
            testCase.verifyEqual(r.maxTopExcessHead_m,zeros(height(r),1),'AbsTol',1e-14);
            testCase.verifyGreaterThan(r.meanComputedTopSettlement_m,zeros(height(r),1));
            testCase.verifyEqual(r.meanComputedTopSettlement_m(end),r.exactTopSettlement_m(end), ...
                'RelTol',.02,'AbsTol',1e-9);
        end
    end
end
