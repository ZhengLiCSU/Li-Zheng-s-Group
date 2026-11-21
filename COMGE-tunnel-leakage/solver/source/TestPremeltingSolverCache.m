classdef TestPremeltingSolverCache < matlab.unittest.TestCase
    methods (Test)
        function cachePreservesStateAndExactTangent(tc)
            [m,p,f,z]=TestPremeltingSolverCache.fixture();
            cached=prepare_premelting_solver(m,p); n=numel(m.weight);
            [a,da]=premelting_switch_ring_state(z,m,p,f,zeros(n,1));
            [b,db]=premelting_switch_ring_state(z,cached,p,f,zeros(n,1));
            tc.verifyEqual(b,a); tc.verifyEqual(db,da);
            tc.verifyEqual(prepare_premelting_solver(cached,p),cached);
        end
        function physicalKeyChangeInvalidatesCache(tc)
            [m,p]=TestPremeltingSolverCache.fixture();
            cached=prepare_premelting_solver(m,p);
            cached.thermalLoad=1.01*cached.thermalLoad;
            updated=prepare_premelting_solver(cached,p);
            tc.verifyNotEqual(updated.elasticSolverCache.key,cached.elasticSolverCache.key);
            tc.verifyEqual(updated.elasticSolverCache.temperatureDerivative, ...
                updated.elasticSolverCache.tangent\updated.thermalLoad);
        end
        function unequalContactOrDamageCannotReuseLinearCache(tc)
            [m,p]=TestPremeltingSolverCache.fixture(); m=prepare_premelting_solver(m,p);
            m.cohesive.contactStiffness_Pam=2*m.cohesive.tensileStiffness_Pam;
            changed=prepare_premelting_solver(m,p); tc.verifyFalse(isfield(changed,'elasticSolverCache'));
            m.cohesive.enabled=true;
            changed=prepare_premelting_solver(m,p); tc.verifyFalse(isfield(changed,'elasticSolverCache'));
        end
    end
    methods (Static,Access=private)
        function [m,p,f,z]=fixture()
            p=joint_threephase_parameters(); f=premelting_switch_parameters();
            c=joint_cohesive_parameters(); c.enabled=false;
            full=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            [m,p]=reduce_premelting_cyclic_symmetry(full,p);
            n=numel(m.weight); ns=numel(m.capacity);
            z=[.1*ones(n,1);ones(n,1);zeros(2*n,1);273.15*ones(n,1); ...
                .1*ones(n,1);279.15*ones(ns,1)];
        end
    end
end
