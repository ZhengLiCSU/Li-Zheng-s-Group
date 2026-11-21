classdef TestTHM20YearValidatorFixture < matlab.unittest.TestCase
    %TESTTHM20YEARVALIDATORFIXTURE End-to-end auditor checks on a SMALL FIXTURE.
    % Twenty actual consecutive warm sealed fixture years are computed;
    % this is neither the research mesh nor a cold-region leakage result.
    % All postprocessing variants use independent copies of the fixture
    % archives. No real campaign files or solver/driver code are changed.
    % The complete positive case proves that the auditor can pass genuine
    % full histories; intentionally damaged copies test rejection paths.

    properties (Access=private)
        Evidence
    end

    methods (TestClassSetup)
        function generateAndAuditSmallFixture(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.Evidence=TestTHM20YearValidatorFixture.prepareEvidence(root);
        end
    end

    methods (Test)
        function actualTwentyWarmFixtureYearsPassNumericalAudit(testCase)
            r=testCase.Evidence;
            testCase.verifyTrue(r.complete.passed);
            testCase.verifyEqual(r.complete.completedYears,20);
            testCase.verifyEqual(height(r.complete.yearSummary),20);
            testCase.verifyTrue(all(r.complete.checks.passed));
            testCase.verifyEqual(r.result.lastYear.timeDay,7300,'AbsTol',1e-7);
            testCase.verifyEqual(r.complete.year20.year20AnnualLeakage_L_per_m,0,'AbsTol',0);
            testCase.verifyTrue(r.complete.numericalValidationOnly);
            testCase.verifyFalse(r.complete.fieldValidation);
        end

        function nineteenYearsCannotBecomePassedYearTwenty(testCase)
            r=testCase.Evidence;
            testCase.verifyFalse(r.nineteen.passed);
            testCase.verifyEqual(r.nineteen.completedYears,19);
            testCase.verifyTrue(any(contains(string(r.nineteen.issues),'complete_twenty_years')));
            testCase.verifyEmpty(fieldnames(r.nineteen.year20));
        end

        function namedAggregateTablesAreAccepted(testCase)
            r=testCase.Evidence;
            testCase.verifyTrue(r.tables.passed);
            testCase.verifyTrue(isequaln(r.tables.yearSummary,r.complete.yearSummary));
            testCase.verifyEqual(r.tables.year20.year20ByJoint_L_per_m,zeros(8,1),'AbsTol',0);
        end

        function nonfiniteRawAggregateLedgerIsRejected(testCase)
            r=testCase.Evidence;
            testCase.verifyFalse(r.nonfinite.passed);
            testCase.verifyTrue(any(contains(string(r.nonfinite.issues),'nonfinite')));
        end

        function negativePhysicalElementPermeabilityIsRejected(testCase)
            r=testCase.Evidence;
            testCase.verifyFalse(r.negativeK.passed);
            testCase.verifyTrue(any(contains(string(r.negativeK.issues), ...
                'Negative element permeability prevents a real diffusive-energy factorization.')));
        end

        function foldedQuadIsRejectedButConsistentReverseIsValid(testCase)
            r=testCase.Evidence;
            testCase.verifyFalse(r.folded.passed);
            testCase.verifyTrue(any(contains(string(r.folded.issues),'Q4 thermal element')));
            testCase.verifyTrue(r.reverse.passed);
            testCase.verifyNotEqual(r.reverse.yearSummary.negativeQ4Elements(1), ...
                r.complete.yearSummary.negativeQ4Elements(1));
            testCase.verifyEqual(r.reverse.year20.year20AnnualLeakage_L_per_m,0,'AbsTol',0);
        end
    end

    methods (Static,Access=private)
        function r=prepareEvidence(root)
            % Fresh run is the default. Explicit evidence reuse is solely a
            % development convenience for qualification-only test edits:
            % setenv('THM20_VALIDATOR_FIXTURE_EVIDENCE',absoluteMatPath).
            % Reaudit the original complete fixture and reject reuse after
            % any inspected solver/validator source hash has changed.
            file=getenv('THM20_VALIDATOR_FIXTURE_EVIDENCE');
            if isempty(file)
                r=TestTHM20YearValidatorFixture.exercise(root);return
            end
            assert(isfile(file),'THM20Fixture:Reuse','Explicit fixture evidence MAT does not exist.');
            saved=load(file,'r');
            assert(isfield(saved,'r')&&strcmp(saved.r.status, ...
                'SMALL WARM SEALED FIXTURE; NOT RESEARCH TUNNEL VALIDATION'), ...
                'THM20Fixture:Reuse','Only the labeled small fixture evidence may be reused.');
            r=saved.r;
            current=validate_thm_20year(r.result,struct('makePlots',false,'saveReport',false));
            assert(current.passed&&isequal(current.sourceHash,r.complete.sourceHash), ...
                'THM20Fixture:Reuse','Fixture/source audit changed. Unset the evidence variable and compute a fresh fixture.');
            fprintf('Explicit unchanged-source SMALL FIXTURE evidence reused: %s.\n',file);
        end

        function r=exercise(root)
            parent=fullfile(root,'THM20_VALIDATOR_FIXTURE_checks');
            if ~isfolder(parent),mkdir(parent);end
            folder=tempname(parent);mkdir(folder);
            caseFolder=fullfile(folder,'WARM_SEALED_SMALL_FIXTURE_20YEARS');
            model=struct('mesh',make_thm_verification_mesh(16),'dtInitial',4*86400, ...
                'dtMax',4*86400,'dtMin',86400/4096,'waterTableDepth',-2, ...
                'topSeepage',false,'geothermalFlux',0,'topMin',4,'topMax',4, ...
                'tunnelMin',4,'tunnelMax',4,'groundwaterInflowT',4, ...
                'groundwaterInflowGradient',0,'initial',struct('mode','profile','surfaceT',4,'gradient',0), ...
                'leakage',struct('jointConductance',zeros(1,8)), ...
                'thm',struct('enabled',true),'solver',struct('jacobian','analytic','reuseAcrossSteps',true));
            result=run_TH_segment_year20_thm(struct('outputFolder',caseFolder,'maxYears',20,'model',model));
            options=struct('makePlots',false,'saveReport',false);
            complete=validate_thm_20year(result,options);
            nineteenResult=result;nineteenResult.completedYears=19;
            nineteenResult.coreCheckpoints=result.coreCheckpoints(1:19);
            nineteenResult.annual=result.annual(result.annual(:,1)<=19,:);
            nineteenResult.monthly=result.monthly(result.monthly(:,1)<=19,:);
            nineteenResult.byJoint=result.byJoint(result.byJoint(:,1)<=19,:);
            nineteenResult.daily=result.daily(result.daily(:,1)<=19,:);
            saved=load(fullfile(caseFolder,'year_19_annual_results.mat'),'out');nineteenResult.lastYear=saved.out;
            nineteen=validate_thm_20year(nineteenResult,options);
            tableResult=result;
            tableResult.annual=array2table(result.annual,'VariableNames',result.annualColumns);
            tableResult.monthly=array2table(result.monthly,'VariableNames',result.monthlyColumns);
            tableResult.byJoint=array2table(result.byJoint,'VariableNames',result.jointColumns);
            tableResult.daily=array2table(result.daily,'VariableNames',result.dailyColumns);
            tables=validate_thm_20year(tableResult,options);
            nonfiniteResult=result;nonfiniteResult.annual(1,5)=NaN;
            nonfinite=validate_thm_20year(nonfiniteResult,options);

            negativeResult=TestTHM20YearValidatorFixture.cloneArchives(result,fullfile(folder,'INVALID_NEGATIVE_K_COPY'));
            file=fullfile(negativeResult.config.outputFolder,'year_01_annual_results.mat');
            saved=load(file,'out');out=saved.out;
            out.state.elementK(1)=-abs(out.state.elementK(1));save(file,'out','-v7');
            negativeK=validate_thm_20year(negativeResult,options);
            foldedResult=TestTHM20YearValidatorFixture.cloneArchives(result,fullfile(folder,'INVALID_FOLDED_Q4_COPY'));
            file=fullfile(foldedResult.config.outputFolder,'year_01_annual_results.mat');
            saved=load(file,'out');out=saved.out;
            out.mesh.qTunnelGlobal(1,:)=out.mesh.qTunnelGlobal(1,[1 3 2 4]);save(file,'out','-v7');
            folded=validate_thm_20year(foldedResult,options);
            reverseResult=TestTHM20YearValidatorFixture.cloneArchives(result,fullfile(folder,'VALID_REVERSE_Q4_COPY'));
            file=fullfile(reverseResult.config.outputFolder,'year_01_annual_results.mat');
            saved=load(file,'out');out=saved.out;
            out.mesh.qTunnelGlobal(1,:)=out.mesh.qTunnelGlobal(1,[1 4 3 2]);save(file,'out','-v7');
            reverse=validate_thm_20year(reverseResult,options);
            r=struct('status','SMALL WARM SEALED FIXTURE; NOT RESEARCH TUNNEL VALIDATION', ...
                'outputFolder',folder,'result',result,'complete',complete,'nineteen',nineteen, ...
                'tables',tables,'nonfinite',nonfinite,'negativeK',negativeK,'folded',folded,'reverse',reverse);
            save(fullfile(folder,'FIXTURE_validator_verification.mat'),'r','-v7');
            fprintf('SMALL FIXTURE validator verification folder: %s.\n',folder);
        end

        function copied=cloneArchives(original,destination)
            mkdir(destination);copied=original;copied.config.outputFolder=destination;
            source=original.config.outputFolder;
            for year=1:20
                stem=sprintf('year_%02d',year);
                for suffix={'_annual_results.mat','_checkpoint.mat'}
                    filename=[stem suffix{1}];
                    copyfile(fullfile(source,filename),fullfile(destination,filename));
                end
            end
        end
    end
end
