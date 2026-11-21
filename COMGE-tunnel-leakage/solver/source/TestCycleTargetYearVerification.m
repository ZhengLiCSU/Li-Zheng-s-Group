classdef TestCycleTargetYearVerification < matlab.unittest.TestCase
    % SMALL warm sealed fixture: qualifies zero-volume and archive contracts.
    properties (Access=private)
        Report
        CaseFolder
    end
    methods (TestClassSetup)
        function runSealedCycleEnabledFixture(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            folder=tempname(fullfile(root,'literature_parameters_20261002'));
            testCase.CaseFolder=folder;
            model=struct('mesh',make_thm_verification_mesh(16),'dtInitial',4*86400, ...
                'dtMax',4*86400,'waterTableDepth',-2,'topSeepage',false, ...
                'geothermalFlux',0,'topMin',4,'topMax',4,'tunnelMin',4,'tunnelMax',4, ...
                'groundwaterInflowT',4,'groundwaterInflowGradient',0, ...
                'initial',struct('mode','profile','surfaceT',4,'gradient',0), ...
                'leakage',struct('jointConductance',zeros(1,8), ...
                'jointCycleDamage',struct('enabled',true)), ...
                'solver',struct('jacobian','analytic'));
            run_TH_segment_year20_thm(struct('outputFolder',folder,'maxYears',2,'model',model));
            testCase.Report=verify_thm_cycle_target_year(folder,2);
        end
    end
    methods (Test)
        function zeroVolumeUsesExactGateAndUndefinedRelativeError(testCase)
            r=testCase.Report;
            testCase.verifyTrue(r.passed);
            testCase.verifyTrue(r.bothAnnualLeakageExactlyZero);
            testCase.verifyEqual([r.coarseQ_L_m,r.fineQ_L_m],[0 0],'AbsTol',0);
            testCase.verifyTrue(isnan(r.annualRelativeFineError));
            testCase.verifyTrue(isnan(r.maxMonthlyDifferenceOverFineAnnual));
        end
        function fullInternalStateAndPhysicsAreInherited(testCase)
            r=testCase.Report;
            testCase.verifyTrue(r.initialXExact&&r.initialMemoryExact ...
                &&r.installedReferenceExact&&r.physicalSignatureExact&&r.localCountsExact);
            testCase.verifyTrue(r.cycleAudit.passed);
            testCase.verifyEqual(r.targetYear,2);
        end
        function existingConditionalEvidenceCannotBeOverwritten(testCase)
            testCase.verifyError(@()verify_thm_cycle_target_year(testCase.CaseFolder,2), ...
                'JointCycleRefinement:Exists');
        end
    end
end
