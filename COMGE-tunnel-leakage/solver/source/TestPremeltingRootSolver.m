classdef TestPremeltingRootSolver < matlab.unittest.TestCase
    methods (Test)
        function stableRootsAgreeWithIndependentScalarReference(tc)
            f=premelting_switch_parameters(); lo=f.minimumFilm_m; hi=f.spinodalFilm_m;
            samples=linspace(lo+1e-12,hi-1e-12,301)';
            potential=premelting_film_potential(samples,f); target=potential.pressure_Pa;
            [roots,info]=solve_premelting_film_roots(target,lo,hi,f);
            reference=zeros(size(target));
            for k=1:numel(target)
                reference(k)=fzero(@(h)TestPremeltingRootSolver.referenceBalance(h,target(k),f), ...
                    [lo,hi],optimset('TolX',1e-22));
            end
            tc.verifyLessThan(max(abs(roots-reference)),2e-22);
            tc.verifyLessThan(max(abs(info.pressureDefect_Pa)),1e-3);
            tc.verifyTrue(all(roots>lo & roots<hi));
        end
        function repulsiveExponentialHasAnalyticInverse(tc)
            f=premelting_switch_parameters(); f.hamaker_J=0;
            expected=linspace(1.1e-9,12e-9,96)';
            target=f.inverseLength_minv*f.shortRangeEnergy_Jm2*exp(-f.inverseLength_minv*expected);
            root=solve_premelting_film_roots(target,1e-9,13e-9,f);
            tc.verifyEqual(root,expected,'AbsTol',1e-22);
        end
        function invalidInteriorCandidateIsRejected(tc)
            f=premelting_switch_parameters();
            tc.verifyError(@()solve_premelting_film_roots(1e20,f.minimumFilm_m,f.spinodalFilm_m,f), ...
                'Premelting:RootBracket');
        end
    end
    methods (Static,Access=private)
        function value=referenceBalance(h,target,f)
            potential=premelting_film_potential(h,f); value=potential.pressure_Pa-target;
        end
    end
end
