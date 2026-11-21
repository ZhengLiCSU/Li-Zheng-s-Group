classdef TestPremelting < matlab.unittest.TestCase
    methods (Test)
        function publishedCoefficientsHaveSIUnits(tc)
            f=premelting_film_parameters(4);
            tc.verifyEqual([f.shortRangeEnergy_Jm2,f.inverseLength_minv,f.hamaker_J], ...
                [0.040969,1.66e9,7.1e-23],'RelTol',1e-14);
        end
        function potentialDerivativeGivesDisjoiningPressure(tc)
            f=premelting_film_parameters(); h=3e-9; epsH=1e-13;
            a=premelting_film_potential(h+epsH,f); b=premelting_film_potential(h-epsH,f);
            v=premelting_film_potential(h,f);
            tc.verifyEqual(-(a.energy_Jm2-b.energy_Jm2)/(2*epsH),v.pressure_Pa,'RelTol',1e-7);
            tc.verifyEqual((a.pressure_Pa-b.pressure_Pa)/(2*epsH),v.pressureDerivative_Pam,'RelTol',1e-7);
        end
        function positiveHamakerTailIsAttractive(tc)
            v=premelting_film_potential(1e-7,premelting_film_parameters());
            tc.verifyLessThan(v.pressure_Pa,0);
            tc.verifyLessThan(v.energy_Jm2,0);
        end
        function pressureSplitObeysChemicalEquilibrium(tc)
            e=premelting_phase_equilibrium(2e5,272.15);
            tc.verifyLessThan(abs(e.affinity_Jkg),1e-9);
            tc.verifyEqual(e.pressureSplit_Pa,e.potential.pressure_Pa,'RelTol',1e-11);
            tc.verifyGreaterThan(e.potential.curvature_Jm4,0);
        end
        function generalizedClapeyronHasTheRightSignAndSlope(tc)
            p=joint_threephase_parameters(); e=premelting_phase_equilibrium(0,273.15,p);
            slope=-p.rhoIce*p.latentHeat_Jkg/p.meltingTemperature_K;
            tc.verifyEqual(e.icePressureDerivative(2),slope,'RelTol',2e-4);
            tc.verifyGreaterThan(e.icePressureDerivative(1),0);
        end
        function coolerEquilibriumHasThinnerFilmAndHigherIcePressure(tc)
            e=premelting_phase_equilibrium([2e5;2e5],[272.8;270.15]);
            tc.verifyGreaterThan(e.film_m(1),e.film_m(2));
            tc.verifyLessThan(e.icePressure_Pa(1),e.icePressure_Pa(2));
        end
        function meltAndMolecularLimitsAreRejected(tc)
            tc.verifyError(@()premelting_phase_equilibrium(2e5,279.15),'Premelting:IceAbsent');
            tc.verifyError(@()premelting_phase_equilibrium(2e5,250),'Premelting:MolecularLimit');
        end
        function cavitatingLiquidIsNotClamped(tc)
            tc.verifyError(@()premelting_phase_equilibrium(-101300,272.15),'Premelting:Cavitation');
        end
        function overlappingFilmsAreNotPatchedToAFluidSlot(tc)
            tc.verifyError(@()premelting_slit_state(2e5,272.15,1e-8),'Premelting:IceTopology');
        end
        function twoPhaseVolumesExactlyFillTheSlit(tc)
            s=premelting_slit_state(2e5,272.15,5e-6);
            tc.verifyEqual(sum(s.phaseVolume_m),s.aperture_m,'AbsTol',1e-20);
            tc.verifyGreaterThan(s.phaseMass_kgm2,[0,0]);
            tc.verifyEqual(s.internalEnergy_Jm2,s.freeEnergy_Jm2+272.15*s.entropy_Jm2K,'AbsTol',1e-10);
        end
        function generalizedGibbsIdentityIncludesFilmEnergy(tc)
            [s,d]=premelting_slit_state(2e5,272.15,5e-6);
            rhs=-s.entropy_Jm2K*[0,1,0]-s.equilibrium.icePressure_Pa*[0,0,1]+ ...
                s.equilibrium.chemicalPotential_Jkg*d.mass;
            tc.verifyEqual(d.freeEnergy,rhs,'RelTol',1e-11,'AbsTol',1e-12);
        end
        function constitutiveJacobianMatchesIndependentDifferences(tc)
            [fd,analytic]=TestPremelting.constitutiveDifferences();
            tc.verifyEqual(fd,analytic,'RelTol',2e-5,'AbsTol',1e-7);
        end
        function stationaryIceSplitsFlowIntoTwoFilms(tc)
            p=joint_threephase_parameters(); s=premelting_slit_state(2e5,272.15,5e-6,p);
            erroneous=s.phaseVolume_m(1)^3/(12*p.viscosity_Pas);
            tc.verifyEqual(erroneous/s.volumeConductance_m3Pas,4,'RelTol',1e-14);
        end
        function finiteElementWallLoadsUseIcePressure(tc)
            [m,p,pl,T,Ts,h]=TestPremelting.ringFixture();
            s=premelting_ring_state(pl,T,Ts,m,p,premelting_film_parameters(),h);
            wrong=solve_cohesive_ring_mechanics(pl,Ts,m,p,h);
            tc.verifyLessThan(norm(s.mechanicalResidual_Nm,inf),1e-4);
            tc.verifyGreaterThan(norm(s.u-wrong.u),1e-7);
        end
        function fullQ4ConstitutiveBridgeHasAConsistentTangent(tc)
            [fd,analytic]=TestPremelting.ringDifference();
            tc.verifyEqual(fd,analytic,'RelTol',2e-5,'AbsTol',1e-8);
        end
        function filmHydraulicPortRetainsPositiveEntropy(tc)
            p=joint_threephase_parameters(); s=premelting_slit_state(2e5,272.15,5e-6,p);
            [~,~,a]=joint_nonisothermal_face([2e5,272.15],[1e5,273],s.chemicalMobility/0.1,1,p);
            tc.verifyGreaterThan(a.entropyProduction_WmK,0);
            tc.verifyEqual(a.entropyProduction_WmK,a.directEntropyProduction_WmK,'RelTol',1e-10);
        end
        function prescribedCoolingConservesEveryClosedCellInventory(tc)
            [a,old,s]=TestPremelting.coolingStep(true);
            tc.verifyLessThan(max(abs(a.localMassDefect_kgm)),1e-12);
            tc.verifyEqual(s.slit.mass_kgm2,old.slit.mass_kgm2,'AbsTol',1e-13);
            tc.verifyLessThan(abs(a.mechanicalWorkQuadratureDefect_Jm),1e-9);
        end
        function limitedSupplyDoesNotImposeTheOpenReservoirIcePressure(tc)
            [a,old,s]=TestPremelting.coolingStep(false);
            tc.verifyLessThan(max(abs(a.localMassDefect_kgm)),1e-12);
            tc.verifyLessThan(max(s.slit.equilibrium.icePressure_Pa-old.slit.equilibrium.icePressure_Pa),20);
            tc.verifyLessThan(max(s.slit.equilibrium.liquidPressure_Pa),1.5e5);
        end
    end
    methods (Static)
        function [a,old,s]=coolingStep(sealed)
            p=joint_threephase_parameters(); f=premelting_film_parameters();
            c=joint_cohesive_parameters(struct('enabled',false));
            m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            n=numel(m.weight); T=272.15*ones(n,1); Ts=279.15*ones(numel(m.capacity),1);
            old=premelting_ring_state(2e5*ones(n,1),T,Ts,m,p,f,zeros(n,1));
            bc=thermo_interface_boundary(m,279.15,279.15,2e5,sealed);
            [s,a]=premelting_prescribed_temperature_step(old,T-0.05,86400,m,bc,p,f);
        end
        function [fd,analytic]=constitutiveDifferences()
            x=[2e5,272.15,5e-6]; delta=[10,1e-5,1e-10]; fd=zeros(3,3);
            [~,d]=premelting_slit_state(x(1),x(2),x(3));
            analytic=[d.mass;d.energy;d.freeEnergy];
            for j=1:3
                xp=x; xm=x; xp(j)=xp(j)+delta(j); xm(j)=xm(j)-delta(j);
                a=premelting_slit_state(xp(1),xp(2),xp(3)); b=premelting_slit_state(xm(1),xm(2),xm(3));
                fd(:,j)=([a.mass_kgm2;a.internalEnergy_Jm2;a.freeEnergy_Jm2]- ...
                    [b.mass_kgm2;b.internalEnergy_Jm2;b.freeEnergy_Jm2])/(2*delta(j));
            end
        end
        function [m,p,pl,T,Ts,h]=ringFixture()
            p=joint_threephase_parameters(); c=joint_cohesive_parameters(struct('strength_Pa',3e5));
            m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            n=numel(m.weight); pl=2e5*ones(n,1); T=272.15*ones(n,1);
            Ts=m.referenceTemperature_K*ones(numel(m.capacity),1); h=zeros(n,1);
        end
        function [fd,analytic]=ringDifference()
            [m,p,pl,T,Ts,h]=TestPremelting.ringFixture(); n=numel(pl); ns=numel(Ts);
            f=premelting_film_parameters();
            [~,d]=premelting_ring_state(pl,T,Ts,m,p,f,h);
            direction=[1e4*cos((1:n)');0.05*sin((1:n)');0.03*cos((1:ns)')];
            x=[pl;T;Ts]; delta=1e-4;
            a=TestPremelting.ringVector(x+delta*direction,m,p,f,h);
            b=TestPremelting.ringVector(x-delta*direction,m,p,f,h);
            fd=(a-b)/(2*delta); analytic=[d.u;d.mass;d.energy]*direction;
        end
        function y=ringVector(x,m,p,f,h)
            n=numel(h); s=premelting_ring_state(x(1:n),x(n+(1:n)),x(2*n+1:end),m,p,f,h);
            y=[s.u;s.slit.mass_kgm2;s.slit.internalEnergy_Jm2];
        end
    end
end
