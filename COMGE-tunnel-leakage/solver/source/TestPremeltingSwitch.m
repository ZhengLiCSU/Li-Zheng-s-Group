classdef TestPremeltingSwitch < matlab.unittest.TestCase
    % Equilibrium thermodynamics, transitions and FE coupling, not calibration.
    methods (Test)
        function phasesShareChemicalPotentialButNotPressure(tc)
            q=premelting_switch_eos([2e5;-2e5],[271;279]);
            tc.verifyLessThan(max(abs(q.g-q.g(:,1)),[],'all'),1e-8);
            tc.verifyGreaterThan(q.pressure(1,2),q.pressure(1,1));
            tc.verifyGreaterThan(q.pressure(:,3)+101325,zeros(2,1));
        end
        function liquidTensionIsAllowedWithoutNegativeVaporPressure(tc)
            q=premelting_switch_eos(-5e5,279.15);
            tc.verifyLessThan(q.pressure(1)+101325,0);
            tc.verifyGreaterThan(q.pressure(3)+101325,0);
        end
        function sameBulkPotentialsAsExistingBranch(tc)
            p=joint_threephase_parameters(); pl=2e5; T=272;
            q=premelting_switch_eos(pl,T,p);
            liquid=joint_threephase_gibbs(pl,T,p);
            ice=joint_threephase_gibbs(q.pressure(2),T,p);
            tc.verifyEqual(q.g(1),liquid.g(1),'AbsTol',1e-10);
            tc.verifyEqual(q.energy(2),ice.energy(2),'AbsTol',1e-9);
        end
        function phaseVolumesAndEnergyClose(tc)
            s=premelting_switch_state(2e5,272.15,5e-6,[.2,.5,.3]);
            tc.verifyEqual(sum(s.phaseVolume_m),5e-6,'AbsTol',1e-20);
            tc.verifyEqual(s.internalEnergy_Jm2,s.freeEnergy_Jm2+272.15*s.entropy_Jm2K,'AbsTol',1e-10);
            tc.verifyGreaterThan(s.phaseMass_kgm2,[0 0 0]);
        end
        function fullGibbsIdentityIncludesAreaWork(tc)
            [s,d]=premelting_switch_state(2e5,272.15,5e-6,[.2,.5,.3]);
            rhs=s.eos.chemicalPotential_Jkg*d.mass-s.entropy_Jm2K*[0 1 0 0 0 0]- ...
                s.wallPressure_Pa*[0 0 1 0 0 0];
            rhs(4:6)=rhs(4:6)+s.grandPotential_Jm2;
            tc.verifyEqual(d.freeEnergy,rhs,'AbsTol',1e-7,'RelTol',1e-10);
        end
        function liquidLimitHasNoGhostIceEnergy(tc)
            p=joint_threephase_parameters(); s=premelting_switch_state(2e5,279.15,5e-6,[1 0 0]);
            tc.verifyEqual(s.surfaceEnergy_Jm2,0);
            tc.verifyEqual(s.phaseMass_kgm2(2:3),[0 0]);
            tc.verifyEqual(s.volumeConductance_m3Pas,(5e-6)^3/(12*p.viscosity_Pas),'RelTol',1e-14);
            tc.verifyEqual(s.wallPressure_Pa,2e5,'AbsTol',1e-10);
        end
        function pureVaporDoesNotNeedCondensateFloor(tc)
            s=premelting_switch_state(-5e5,279.15,5e-6,[0 0 1]);
            tc.verifyEqual(s.phaseMass_kgm2(1:2),[0 0]);
            tc.verifyGreaterThan(s.phaseMass_kgm2(3),0);
            tc.verifyEqual(s.chemicalMobility,0);
            tc.verifyEqual(phase_simplex_residual([0 0 1],s.grandPotential_Jm2,6),[0 0 0],'AbsTol',1e-12);
        end
        function wetIceLimitRetainsTwoSeparateFilms(tc)
            p=joint_threephase_parameters(); f=premelting_switch_parameters();
            a=premelting_switch_state(2e5,272.15,5e-6,[0 1 0],p,f);
            b=premelting_slit_state(2e5,272.15,5e-6,p,f);
            tc.verifyEqual(a.mass_kgm2,b.mass_kgm2,'RelTol',1e-13);
            tc.verifyEqual(a.internalEnergy_Jm2,b.internalEnergy_Jm2+2*f.iceWaterSurfaceEnergy_Jm2,'AbsTol',1e-10);
            tc.verifyEqual(a.chemicalMobility,b.chemicalMobility,'RelTol',1e-12);
        end
        function allThreePureStatesHaveCorrectKktSign(tc)
            warm=premelting_switch_state(2e5,279.15,5e-6,[1 0 0]);
            cold=premelting_switch_state(2e5,272.15,5e-6,[0 1 0]);
            tc.verifyEqual(phase_simplex_residual([1 0 0],warm.grandPotential_Jm2,6),[0 0 0],'AbsTol',1e-12);
            tc.verifyEqual(phase_simplex_residual([0 1 0],cold.grandPotential_Jm2,6),[0 0 0],'AbsTol',1e-12);
            tc.verifyGreaterThan(norm(phase_simplex_residual([1 0 0],cold.grandPotential_Jm2,6)),0.1);
        end
        function simplexResidualRejectsNonunitArea(tc)
            r=phase_simplex_residual([.2 .3 .4],[0 0 0],6);
            tc.verifyEqual(sum(r),-.1,'AbsTol',1e-14);
        end
        function interfaceEnergyVanishesContinuouslyWithIceArea(tc)
            a=premelting_switch_state(2e5,272.15,5e-6,[1-1e-5 1e-5 0]);
            b=premelting_switch_state(2e5,272.15,5e-6,[1-2e-5 2e-5 0]);
            tc.verifyEqual(b.surfaceEnergy_Jm2,2*a.surfaceEnergy_Jm2,'RelTol',1e-14);
        end
        function exactConstitutiveDerivativesMatchDifferences(tc)
            error=TestPremeltingSwitch.localDerivativeError();
            tc.verifyLessThan(error,2e-6);
        end
        function monolithicJacobianIncludesWallAndPhaseCoupling(tc)
            error=TestPremeltingSwitch.ringDerivativeError();
            tc.verifyLessThan(error,2e-6);
        end
        function sealedEquilibriumDoesNotCreateWaterOrHeat(tc)
            [m,p,f,z]=TestPremeltingSwitch.fixture(); n=numel(m.weight);
            bc=thermo_interface_boundary(m,279.15,279.15,1e5,true);
            [next,a]=premelting_switch_ring_step(z,86400,m,bc,p,f,zeros(n,1));
            tc.verifyEqual(next,z,'AbsTol',1e-10);
            tc.verifyLessThan(max(abs(a.localMassDefect_kgm)),1e-12);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
        end
        function pureVaporHeatingNeedsNoArtificialLiquidInventory(tc)
            [m,p,f,z]=TestPremeltingSwitch.fixture(); n=numel(m.weight);
            eos=premelting_switch_eos(-5e5,279.15,p);
            z(1:n)=-.5; z(n+(1:2*n))=0; z(3*n+(1:n))=1;
            z(5*n+(1:n))=eos.pressure(3)/1e6;
            bc=thermo_interface_boundary(m,280.15,280.15,1e5,true);
            old=premelting_switch_ring_state(z,m,p,f,zeros(n,1));
            [~,a]=premelting_switch_ring_step(z,3600,m,bc,p,f,zeros(n,1));
            tc.verifyEqual(a.state.areaFraction,repmat([0 0 1],n,1),'AbsTol',1e-12);
            tc.verifyEqual(sum(m.weight.*a.state.mass_kgm2),sum(m.weight.*old.mass_kgm2),'AbsTol',1e-12);
            tc.verifyGreaterThan(mean(a.state.x(:,4)),279.15);
            tc.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
        end
    end
    methods (Static, Access=private)
        function e=localDerivativeError()
            x=[2e5,272.15,5e-6,.2,.5,.3]; step=[.1,1e-5,1e-11,1e-6,1e-6,1e-6];
            [~,d]=premelting_switch_state(x(1),x(2),x(3),x(4:6));
            analytic=[d.mass;d.energy;d.entropy;d.wallPressure;d.chemicalMobility];
            fd=zeros(5,6);
            for k=1:6
                a=x; b=x; a(k)=a(k)+step(k); b(k)=b(k)-step(k);
                sa=premelting_switch_state(a(1),a(2),a(3),a(4:6));
                sb=premelting_switch_state(b(1),b(2),b(3),b(4:6));
                fd(:,k)=(TestPremeltingSwitch.values(sa)-TestPremeltingSwitch.values(sb))/(2*step(k));
            end
            e=max(abs(fd-analytic)./max(abs(analytic),max(abs(analytic),[],2)*1e-8),[],'all');
        end
        function v=values(s)
            v=[s.mass_kgm2;s.internalEnergy_Jm2;s.entropy_Jm2K;s.wallPressure_Pa;s.chemicalMobility];
        end
        function e=ringDerivativeError()
            [m,p,f,z]=TestPremeltingSwitch.fixture(); n=numel(m.weight); x=z;
            x(n+(1:3*n))=repmat([.2 .5 .3],n,1); x(4*n+(1:n))=272.15;
            bc=thermo_interface_boundary(m,279.15,269.15,2e5,false);
            [~,J]=premelting_switch_ring_residual(x,z,100,m,bc,p,f,zeros(n,1));
            direction=sin((1:numel(x))'); direction(6*n+1:end)=.01*direction(6*n+1:end); h=1e-5;
            rp=premelting_switch_ring_residual(x+h*direction,z,100,m,bc,p,f,zeros(n,1));
            rm=premelting_switch_ring_residual(x-h*direction,z,100,m,bc,p,f,zeros(n,1));
            fd=(rp-rm)/(2*h); e=norm(fd-J*direction)/norm(fd);
        end
        function [m,p,f,z]=fixture()
            p=joint_threephase_parameters(); f=premelting_switch_parameters(); c=joint_cohesive_parameters(); c.enabled=false;
            m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            n=numel(m.weight); ns=numel(m.capacity);
            z=[.1*ones(n,1);ones(n,1);zeros(2*n,1);279.15*ones(n,1);.1*ones(n,1);279.15*ones(ns,1)];
        end
    end
end
