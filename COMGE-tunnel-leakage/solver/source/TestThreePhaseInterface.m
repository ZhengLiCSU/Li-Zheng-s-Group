classdef TestThreePhaseInterface < matlab.unittest.TestCase
    methods (Test)
        function allChemicalPotentialsMeetAtTriplePoint(tc)
            p=joint_threephase_parameters();
            q=joint_threephase_gibbs(p.triplePressure_Pa-p.atmosphericPressure_Pa,273.16,p);
            tc.verifyLessThan(max(q.g)-min(q.g),1e-7);
        end
        function sublimationMatchesIndependentIAPWSWithinOnePercent(tc)
            T=(250:2:272)'; predicted=arrayfun(@threephase_saturation_pressure,T);
            reference=iapws_sublimation_reference(T);
            tc.verifyLessThan(max(abs(predicted./reference-1)),0.01);
        end
        function gibbsDerivativesGiveVolumesAndEntropies(tc)
            p=joint_threephase_parameters(); pg=-1e5; T=270; h=0.01;
            q=joint_threephase_gibbs(pg,T,p);
            a=joint_threephase_gibbs(pg+h,T,p); b=joint_threephase_gibbs(pg-h,T,p);
            c=joint_threephase_gibbs(pg,T+h,p); d=joint_threephase_gibbs(pg,T-h,p);
            tc.verifyEqual((a.g-b.g)/(2*h),q.volume,'RelTol',2e-5);
            tc.verifyEqual(-(c.g-d.g)/(2*h),q.entropy,'RelTol',1e-6);
        end
        function energyDerivativesIncludeAmbientWork(tc)
            p=joint_threephase_parameters(); pg=-1e5; T=270; h=0.01;
            q=joint_threephase_gibbs(pg,T,p);
            a=joint_threephase_gibbs(pg+h,T,p); b=joint_threephase_gibbs(pg-h,T,p);
            c=joint_threephase_gibbs(pg,T+h,p); d=joint_threephase_gibbs(pg,T-h,p);
            tc.verifyEqual((a.energy-b.energy)/(2*h),q.energyPressure,'RelTol',3e-5,'AbsTol',1e-8);
            tc.verifyEqual((c.energy-d.energy)/(2*h),q.energyTemperature,'RelTol',1e-6);
        end
        function zeroAbsolutePressureIsRejected(tc)
            p=joint_threephase_parameters();
            tc.verifyError(@()joint_threephase_gibbs(-p.atmosphericPressure_Pa,270,p),'ThreePhase:Domain');
        end
        function stateVolumeClosesWithNonzeroVaporMass(tc)
            [m,p,z]=TestThreePhaseInterface.fixture(); n=numel(m.weight);
            z(1:n)=log(500/p.atmosphericPressure_Pa); z(n+(1:n))=1; z(2*n+(1:n))=0.15;
            s=threephase_interface_state(z,m,p);
            tc.verifyEqual(s.liquidVolume_m+s.iceVolume_m+s.vaporVolume_m,s.aperture_m,'RelTol',1e-13);
            tc.verifyGreaterThan(s.phaseMass_kgm2(:,3),zeros(n,1));
            tc.verifyEqual(sum(s.phaseMass_kgm2,2),s.x(:,2));
        end
        function fullTangentMatchesIndependentDirectionalDifference(tc)
            [m,p,z]=TestThreePhaseInterface.fixture(); n=numel(m.weight);
            z(1:n)=log(500/p.atmosphericPressure_Pa); z(n+(1:n))=0.4; z(2*n+(1:n))=0.15; z(3*n+(1:n))=270;
            bc=thermo_interface_boundary(m,276,274);
            [~,J]=threephase_interface_residual(z,z,60,m,bc,p);
            direction=[0.1*cos((1:n)');0.01*sin((1:2*n)');0.1*cos((1:numel(z)-3*n)')];
            h=1e-5;
            fd=(threephase_interface_residual(z+h*direction,z,60,m,bc,p)- ...
                threephase_interface_residual(z-h*direction,z,60,m,bc,p))/(2*h);
            tc.verifyEqual(J*direction,fd,'RelTol',4e-5,'AbsTol',3e-6);
        end
        function warmSaturatedFluxAgreesWithExistingBranch(tc)
            [m,p,z]=TestThreePhaseInterface.fixture(); n=numel(m.weight);
            old=[p.atmosphericPressure_Pa*expm1(z(1:n));z(n+(1:n));z(3*n+1:end)];
            bc=thermo_interface_boundary(m,279.15,279.15);
            [zn,a]=threephase_interface_step(z,60,m,bc,p);
            [ref,b]=thermo_interface_step(old,60,m,bc,p);
            tc.verifyEqual(zn(2*n+(1:n)),zeros(n,1),'AbsTol',1e-14);
            tc.verifyEqual(a.faceFlux,b.faceFlux,'RelTol',1e-5,'AbsTol',1e-10);
            tc.verifyEqual(a.state.pressure,ref(1:n),'AbsTol',1e-3);
        end
        function closedWarmVaporCondensesConservatively(tc)
            [m,p,z]=TestThreePhaseInterface.fixture(); n=numel(m.weight);
            z(2*n+(1:n))=0.1; m.bathConductance(:)=0;
            before=threephase_interface_state(z,m,p);
            bc=thermo_interface_boundary(m,279.15,279.15,0,true); bc.insulatedHeat(:)=true;
            [zn,a]=threephase_interface_step(z,60,m,bc,p);
            tc.verifyLessThan(max(abs(a.localMassDefect_kgm)),1e-10);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
            tc.verifyLessThan(max(zn(2*n+(1:n))),0.1);
            tc.verifyEqual(sum(m.weight.*a.state.x(:,2)),sum(m.weight.*before.x(:,2)),'AbsTol',1e-11);
            tc.verifyGreaterThan(a.entropyProduced_JmK,0);
        end
        function frozenVaporEquilibriumIsStationary(tc)
            [m,p,z]=TestThreePhaseInterface.fixture(); n=numel(m.weight);
            T=270; pv=threephase_saturation_pressure(T,p);
            z(1:n)=log(pv/p.atmosphericPressure_Pa); z(n+(1:n))=1;
            z(2*n+(1:n))=0.15; z(3*n+1:end)=T;
            bc=thermo_interface_boundary(m,T,T,0,true); bc.insulatedHeat(:)=true;
            [zn,a]=threephase_interface_step(z,60,m,bc,p);
            tc.verifyEqual(zn,z,'AbsTol',1e-8);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),1e-8);
        end
        function steadyHeatLedgerCannotBeSkippedByScaledTolerance(tc)
            % Reproduce the late-time near-steady regime without a saved
            % fixture or weaker tolerance: long warm run reaches that state.
            [m,p,~]=TestThreePhaseInterface.fixture();
            out=solve_threephase_ring(m,p,struct('cycles',0.6,'temperatureAmplitude_K',0, ...
                'dtDays',0.0625,'maximumPhaseIncrement',0.0125, ...
                'maximumPressureIncrement_Pa',5000,'maximumLogPressureIncrement',0.125,'verbose',false));
            tc.verifyTrue(out.completed);
            tc.verifyLessThan(abs(out.audit(:,4)),out.audit(:,8));
        end
    end
    methods (Static, Access=private)
        function [m,p,z]=fixture()
            p=joint_threephase_parameters();
            m=build_thermo_interface_ring(make_thm_verification_mesh(16),3,p);
            n=numel(m.weight); ns=numel(m.capacity);
            z=[log1p(1e5/p.atmosphericPressure_Pa)*ones(n,1);zeros(2*n,1);279.15*ones(n+ns,1)];
        end
    end
end
