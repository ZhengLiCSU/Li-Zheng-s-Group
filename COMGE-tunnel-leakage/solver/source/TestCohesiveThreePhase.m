classdef TestCohesiveThreePhase < matlab.unittest.TestCase
    methods (Test)
        function completeSeparationCostsSpecifiedFractureEnergy(tc)
            c=joint_cohesive_parameters();
            s=joint_cohesive_energy(2*c.failureGap_m,0,c);
            tc.verifyEqual(s.surfaceEnergy_Jm2,c.fractureEnergy_Jm2);
            tc.verifyEqual(s.elasticEnergy_Jm2,0);
            tc.verifyEqual(s.traction_Pa,0);
            tc.verifyEqual(s.damage,1);
        end
        function compressiveContactSurvivesCompleteTensileFailure(tc)
            c=joint_cohesive_parameters();
            s=joint_cohesive_energy(-1e-7,2*c.failureGap_m,c);
            tc.verifyEqual(s.traction_Pa,-1e-7*c.contactStiffness_Pam);
            tc.verifyEqual(s.tangent_Pam,c.contactStiffness_Pam);
            tc.verifyEqual(s.surfaceEnergy_Jm2,c.fractureEnergy_Jm2);
        end
        function unloadingDoesNotHealOrCreateFurtherDamage(tc)
            c=joint_cohesive_parameters(); peak=0.4*c.failureGap_m;
            loaded=joint_cohesive_energy(peak,0,c);
            unloaded=joint_cohesive_energy(peak/5,loaded.kappa_m,c);
            tc.verifyEqual(unloaded.damage,loaded.damage);
            tc.verifyEqual(unloaded.surfaceEnergy_Jm2,loaded.surfaceEnergy_Jm2);
            tc.verifyEqual(unloaded.traction_Pa,loaded.traction_Pa/5,'RelTol',1e-12);
        end
        function monotonicWorkMatchesAnalyticalTriangleArea(tc)
            c=joint_cohesive_parameters(); g=0.4*c.failureGap_m;
            s=joint_cohesive_energy(g,0,c);
            expected=0.5*c.strength_Pa*c.onsetGap_m+ ...
                0.5*(c.strength_Pa+s.traction_Pa)*(g-c.onsetGap_m);
            tc.verifyEqual(s.elasticEnergy_Jm2+s.surfaceEnergy_Jm2,expected,'RelTol',1e-13);
        end
        function fractureDrivingForceMatchesSurfaceEnergyDerivative(tc)
            c=joint_cohesive_parameters(); g=0.4*c.failureGap_m;
            s=joint_cohesive_energy(g,0,c);
            tc.verifyEqual(s.surfaceDerivative/s.damageDerivative,0.5*c.tensileStiffness_Pam*g^2,'RelTol',1e-13);
        end
        function reproducesExistingPublishedBilinearModeOneLimit(tc)
            c=joint_cohesive_parameters(); g=0.4*c.failureGap_m;
            material=struct('strength_Pa',c.strength_Pa,'crackingGap_m',c.onsetGap_m,'debondingGap_m',c.failureGap_m);
            [traction,tangent]=luo2026_cohesive_mode1(g,0,struct(),material);
            s=joint_cohesive_energy(g,0,c);
            tc.verifyEqual(s.traction_Pa,traction,'RelTol',1e-12);
            tc.verifyEqual(s.tangent_Pam,tangent,'RelTol',1e-12);
        end
        function activeSofteningTangentMatchesFiniteDifference(tc)
            c=joint_cohesive_parameters(); g=0.2*c.failureGap_m; h=1e-9;
            s=joint_cohesive_energy(g,0,c);
            a=joint_cohesive_energy(g+h,0,c); b=joint_cohesive_energy(g-h,0,c);
            tc.verifyEqual((a.traction_Pa-b.traction_Pa)/(2*h),s.tangent_Pam,'RelTol',1e-7);
        end
        function unstableSupportedMaterialIsNotArtificiallyStabilized(tc)
            p=joint_threephase_parameters(); c=joint_cohesive_parameters();
            tc.verifyError(@()build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c),'Cohesive:Stability');
        end
        function disabledDamageDoesNotApplyANonexistentSofteningBranch(tc)
            p=joint_threephase_parameters(); c=joint_cohesive_parameters(struct('enabled',false));
            m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            tc.verifyFalse(m.cohesive.enabled);
        end
        function selectedKinkDerivativeAlsoMatchesFractureEnergyMetadata(tc)
            [m,p,z,h]=TestCohesiveThreePhase.fixture(); n=numel(h);
            Ts=277.15+zeros(numel(z)-4*n,1); pg=1e5+zeros(n,1);
            a=solve_cohesive_ring_mechanics(pg,Ts,m,p,h);
            history=a.local.kappa_m;
            loading=solve_cohesive_ring_mechanics(pg,Ts,m,p,history,true(n,1));
            unloading=solve_cohesive_ring_mechanics(pg,Ts,m,p,history,false(n,1));
            tc.verifyEqual(loading.u,unloading.u,'AbsTol',1e-15);
            tc.verifyGreaterThan(loading.local.damageDerivative,zeros(n,1));
            tc.verifyEqual(unloading.local.damageDerivative,zeros(n,1));
            tc.verifyEqual(loading.local.surfaceDerivative./loading.local.damageDerivative, ...
                .5*m.cohesive.tensileStiffness_Pam*history.^2,'RelTol',1e-12);
        end
        function undamagedHostRecoversPreviousRing(tc)
            [m,p,z,h]=TestCohesiveThreePhase.fixture();
            old=build_thermo_interface_ring(make_thm_verification_mesh(16),3,p);
            s=cohesive_threephase_state(z,m,p,h); ref=threephase_interface_state(z,old,p);
            tc.verifyEqual(s.u,ref.u,'AbsTol',1e-15);
            tc.verifyEqual(s.x,ref.x,'RelTol',1e-10,'AbsTol',1e-14);
        end
        function nonlinearHostHasConsistentPressureTangent(tc)
            [m,p,z,h]=TestCohesiveThreePhase.fixture(); n=numel(m.weight); Ts=z(4*n+1:end)-2;
            pg=1e5*ones(n,1); direction=1e5*cos((1:n)');
            a=solve_cohesive_ring_mechanics(pg,Ts,m,p,h);
            plus=solve_cohesive_ring_mechanics(pg+1e-5*direction,Ts,m,p,h);
            minus=solve_cohesive_ring_mechanics(pg-1e-5*direction,Ts,m,p,h);
            tc.verifyEqual((plus.u-minus.u)/2e-5,a.pressureDerivative*direction,'RelTol',1e-5,'AbsTol',1e-12);
        end
        function fullCoupledTangentSurvivesDamageAndThreePhases(tc)
            [m,p,z,h]=TestCohesiveThreePhase.fixture(); n=numel(m.weight);
            z(n+(1:n))=0.4; z(2*n+(1:n))=0.1; z(3*n+(1:n))=273.1; z(4*n+1:end)=277.15;
            bc=thermo_interface_boundary(m,279.15,275);
            [~,J]=cohesive_threephase_residual(z,z,60,m,bc,p,h);
            direction=[0.1*cos((1:n)');0.01*sin((1:2*n)');0.1*cos((1:numel(z)-3*n)')];
            step=1e-5;
            fd=(cohesive_threephase_residual(z+step*direction,z,60,m,bc,p,h)- ...
                cohesive_threephase_residual(z-step*direction,z,60,m,bc,p,h))/(2*step);
            tc.verifyEqual(J*direction,fd,'RelTol',5e-5,'AbsTol',1e-5);
        end
        function priorDamageIncreasesLoadedOpeningWithoutLeakMultiplier(tc)
            [m,p,z,h]=TestCohesiveThreePhase.fixture();
            intact=cohesive_threephase_state(z,m,p,h);
            damaged=cohesive_threephase_state(z,m,p,1e-5+zeros(size(h)));
            tc.verifyGreaterThan(damaged.aperture_m,intact.aperture_m);
            tc.verifyGreaterThan(damaged.liquidVolume_m.^3,intact.liquidVolume_m.^3);
            tc.verifyEqual(damaged.pressure,intact.pressure);
        end
        function acceptedWarmStepConservesMassAndTotalEnergy(tc)
            [m,p,z,h]=TestCohesiveThreePhase.fixture();
            [~,a]=cohesive_threephase_step(z,60,m,thermo_interface_boundary(m,279.15,279.15),p,h);
            tc.verifyLessThan(max(abs(a.localMassDefect_kgm)),1e-10);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
            tc.verifyGreaterThanOrEqual(a.fractureSurfaceEnergyChange_Jm,0);
            tc.verifyLessThan(norm(a.mechanicalResidual_Nm,inf),1e-4);
        end
        function historicalMaximumKinkUsesDirectionalTangent(tc)
            f=load(fullfile(fileparts(mfilename('fullpath')),'cohesive_kink_regression.mat'));
            [next,a]=cohesive_threephase_step(f.old,f.dt,f.model,f.bc,f.p,f.history);
            tc.verifyTrue(all(isfinite(next)));
            tc.verifyLessThan(a.iterations,10);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
            tc.verifyGreaterThanOrEqual(a.state.cohesive.kappa_m,f.history);
        end
        function fullyFrozenStationaryGapDoesNotToggleOnRoundoff(tc)
            f=load(fullfile(fileparts(mfilename('fullpath')),'cohesive_frozen_kink_regression.mat'));
            [next,a]=cohesive_threephase_step(f.old,f.dt,f.model,f.bc,f.p,f.history);
            n=numel(f.history);
            tc.verifyEqual(next(n+(1:n)),ones(n,1),'AbsTol',1e-12);
            tc.verifyLessThan(a.iterations,10);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
        end
        function warmReloadBalancesFlowWithoutCreatingDamage(tc)
            [m,p,~,h]=TestCohesiveThreePhase.fixture();
            a=cohesive_isothermal_leak_probe(m,p,h);
            tc.verifyGreaterThan(a.rate_Lmday,0);
            tc.verifyEqual(a.maxDamage,0);
            tc.verifyLessThan(norm(a.massResidual_kgms,inf),1e-14);
            tc.verifyEqual(a.jointRate_Lmday,a.rate_Lmday/8+zeros(8,1),'RelTol',1e-9);
        end
    end
    methods (Static, Access=private)
        function [m,p,z,h]=fixture()
            p=joint_threephase_parameters(); c=joint_cohesive_parameters(struct('strength_Pa',3e5));
            m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            n=numel(m.weight); ns=numel(m.capacity); h=zeros(n,1);
            z=[log1p(1e5/p.atmosphericPressure_Pa)*ones(n,1);zeros(2*n,1);279.15*ones(n+ns,1)];
        end
    end
end
