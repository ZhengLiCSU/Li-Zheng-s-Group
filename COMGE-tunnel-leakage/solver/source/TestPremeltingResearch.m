classdef TestPremeltingResearch < matlab.unittest.TestCase
    methods (Test)
        function originalMeshIsExactlyPreserved(tc)
            tc.verifyEqual(make_refined_thm_verification_mesh(16,3),make_thm_verification_mesh(16));
        end
        function cyclicMechanicsEqualsEightSegmentMechanics(tc)
            [m,mr,p,pr,f,z]=TestPremeltingResearch.fixture();
            z(4*mr.nq+(1:mr.nq))=272.15;
            ss=premelting_switch_ring_state(z,mr,pr,f,zeros(mr.nq,1));
            sf=premelting_switch_ring_state(lift_premelting_symmetry_state(z,mr),m,p,f,zeros(8*mr.nq,1));
            tc.verifyEqual(sf.u,mr.symmetry.displacementLift*ss.u,'AbsTol',1e-14);
            tc.verifyEqual(sf.mass_kgm2,mr.symmetry.cellLift*ss.mass_kgm2,'AbsTol',1e-15);
        end
        function cyclicCoupledStepEqualsEightSegmentStep(tc)
            [m,mr,p,pr,f,z]=TestPremeltingResearch.fixture();
            bcr=thermo_interface_boundary(mr,279.15,269.15,2e5,false);
            bcf=thermo_interface_boundary(m,279.15,269.15,2e5,false);
            [zr,ar]=premelting_switch_ring_step(z,3600,mr,bcr,pr,f,zeros(mr.nq,1));
            [zf,af]=premelting_switch_ring_step(lift_premelting_symmetry_state(z,mr),3600,m,bcf,p,f,zeros(8*mr.nq,1));
            tc.verifyEqual(zf,lift_premelting_symmetry_state(zr,mr),'AbsTol',1e-8);
            tc.verifyEqual(af.fluidEnergyChange_Jm,8*ar.fluidEnergyChange_Jm,'AbsTol',1e-6);
        end
        function zeroVarianceFieldPreservesResidualAndJacobian(tc)
            [m,~,p,~,f,zr]=TestPremeltingResearch.fixture();
            mr=reduce_premelting_cyclic_symmetry(m,p); z=lift_premelting_symmetry_state(zr,mr);
            pf=apply_premelting_aperture_field(m,p,struct('xi',[0,1],'aperture_m',p.aperture0_m*ones(8,2)));
            bc=thermo_interface_boundary(m,279.15,279.15,2e5,false);
            [r1,J1]=premelting_switch_ring_residual(z,z,100,m,bc,p,f,zeros(numel(m.weight),1));
            [r2,J2]=premelting_switch_ring_residual(z,z,100,m,bc,pf,f,zeros(numel(m.weight),1));
            tc.verifyEqual(r2,r1,'AbsTol',1e-14); tc.verifyEqual(J2,J1,'AbsTol',1e-12);
        end
        function interpolationIsConsistentAcrossMeshRefinement(tc)
            [m,~,p]=TestPremeltingResearch.fixture();
            c=joint_cohesive_parameters(); c.enabled=false;
            fine=build_cohesive_threephase_ring(make_thm_verification_mesh(16),9,p,c);
            field=struct('xi',[0,1],'aperture_m',repmat([4e-6,6e-6],8,1));
            a=apply_premelting_aperture_field(m,p,field); b=apply_premelting_aperture_field(fine,p,field);
            fineValues=reshape(b.aperture0_m,9,8);
            tc.verifyEqual(reshape(a.aperture0_m,3,8),fineValues([2,5,8],:),'AbsTol',1e-18);
        end
        function heterogeneousAperturesChangeInventoryAndTransportTogether(tc)
            [m,mr,p,~,f,zr]=TestPremeltingResearch.fixture();
            field=struct('xi',[0,1],'aperture_m',repmat(linspace(4e-6,6e-6,8)',1,2));
            pf=apply_premelting_aperture_field(m,p,field); z=lift_premelting_symmetry_state(zr,mr);
            s=premelting_switch_ring_state(z,m,pf,f,zeros(numel(m.weight),1));
            tc.verifyGreaterThan(s.mass_kgm2(end),s.mass_kgm2(1));
            tc.verifyGreaterThan(s.volumeConductance_m3Pas(end),s.volumeConductance_m3Pas(1));
            tc.verifyEqual(s.volumeConductance_m3Pas,s.aperture_m.^3/(12*p.viscosity_Pas),'AbsTol',1e-25);
        end
        function asymmetricFieldCannotUseCyclicShortcut(tc)
            [m,~,p]=TestPremeltingResearch.fixture();
            pf=apply_premelting_aperture_field(m,p,struct('xi',[0,1],'aperture_m',repmat(linspace(4e-6,6e-6,8)',1,2)));
            tc.verifyError(@()reduce_premelting_cyclic_symmetry(m,pf),'PhaseSwitch:Symmetry');
        end
        function heterogeneousResidualJacobianMatchesFiniteDifference(tc)
            [m,mr,p,~,f,zr]=TestPremeltingResearch.fixture();
            pf=apply_premelting_aperture_field(m,p,struct('xi',[0,1],'aperture_m',repmat([4e-6,6e-6],8,1)));
            z=lift_premelting_symmetry_state(zr,mr); bc=thermo_interface_boundary(m,279.15,278.15,2e5,false);
            [~,J]=premelting_switch_ring_residual(z,z,100,m,bc,pf,f,zeros(numel(m.weight),1));
            direction=cos((1:numel(z))'); direction(24+(1:72))=0;
            e=1e-6;
            rp=premelting_switch_ring_residual(z+e*direction,z,100,m,bc,pf,f,zeros(numel(m.weight),1));
            rm=premelting_switch_ring_residual(z-e*direction,z,100,m,bc,pf,f,zeros(numel(m.weight),1));
            tc.verifyLessThan(norm((rp-rm)/(2*e)-J*direction)/norm(J*direction),1e-6);
        end
        function negativeTrialApertureHasRecoverableDomainIdentifier(tc)
            [~,mr,~,pr,f,z]=TestPremeltingResearch.fixture();
            z(5*mr.nq+(1:mr.nq))=-100;
            tc.verifyError(@()premelting_switch_ring_state(z,mr,pr,f,zeros(mr.nq,1)),'PhaseSwitch:Domain');
        end
        function invalidSampleDoesNotBecomeZeroLeakage(tc)
            cfg=premelting_research_case(); cfg.apertureField=struct('xi',[0,1],'aperture_m',-ones(8,2));
            folder=tempname(pwd); root=pwd;
            cleanup=onCleanup(@()TestPremeltingResearch.removeTestFolder(folder,root));
            result=run_premelting_sample(cfg,folder);
            tc.verifyEqual(result.status,'input_or_physical_domain_failure');
            tc.verifyEmpty(result.qoi); tc.verifyTrue(isfile(fullfile(folder,'input.mat')));
            tc.verifyError(@()run_premelting_sample(cfg,folder),'PhaseSwitch:OutputExists');
        end
        function boundarySeriesUsesSIUnitsAndForbidsExtrapolation(tc)
            o=struct('boundarySchedule',struct('time_s',[0,10,20], ...
                'outerTemperature_K',[279,280,281],'innerTemperature_K',[273,270,271], ...
                'reservoirPressure_Pa',[1e5,2e5,3e5]));
            [a,b,h]=premelting_boundary_at(5,279.15,o);
            tc.verifyEqual([a,b,h],[279.5,271.5,1.5e5]);
            tc.verifyError(@()premelting_boundary_at(21,279.15,o),'PhaseSwitch:Input');
        end
        function pressureControlSeesIceAndIgnoresAbsentGhostPressure(tc)
            a=struct('areaFraction',[1,0,0],'eos',struct('pressure',[0,-1e7,0]), ...
                'liquidPressure',0,'pressure',0,'x',[0,1,0,273]);
            b=a; b.eos.pressure(2)=-2e7;
            change=premelting_transition_change(a,b); tc.verifyEqual(change.activeIcePressure_Pa,0);
            b.areaFraction=[.99,.01,0]; change=premelting_transition_change(a,b);
            tc.verifyEqual(change.activeIcePressure_Pa,1e7);
        end
        function solverLandsOnBoundaryKnots(tc)
            [~,mr,~,pr,f]=TestPremeltingResearch.fixture();
            schedule=struct('time_s',[0,10,20],'outerTemperature_K',[279.15,279.15,279.15], ...
                'innerTemperature_K',[279.15,278.15,279.15],'reservoirPressure_Pa',[2e5,2e5,2e5]);
            out=solve_premelting_switch_ring(mr,pr,f,struct('cycles',1,'periodDays',20/86400, ...
                'dtDays',60/86400,'boundarySchedule',schedule,'sealed',true,'verbose',false));
            tc.verifyEqual(out.time_s,[0;10;20]);
        end
        function hydraulicRefinementNeedsResolvedMechanicalTrace(tc)
            [m,~,p]=TestPremeltingResearch.fixture(); evidence=check_premelting_field_mesh(m);
            tc.verifyTrue(evidence.traceFullRowRank);
            c=joint_cohesive_parameters(); c.enabled=false;
            unresolved=build_cohesive_threephase_ring(make_thm_verification_mesh(16),6,p,c);
            tc.verifyError(@()check_premelting_field_mesh(unresolved),'PhaseSwitch:TraceResolution');
            resolved=build_cohesive_threephase_ring(make_refined_thm_verification_mesh(16,6),6,p,c);
            evidence=check_premelting_field_mesh(resolved); tc.verifyTrue(evidence.traceFullRowRank);
        end
        function iapwsTransportJacobianMatchesFiniteDifference(tc)
            [~,mr,~,pr,f,z]=TestPremeltingResearch.fixture(); pr.viscosityModel='iapws2008';
            bc=thermo_interface_boundary(mr,279.15,278.15,2e5,false);
            [~,J]=premelting_switch_ring_residual(z,z,100,mr,bc,pr,f,zeros(mr.nq,1));
            direction=cos((1:numel(z))'); direction(mr.nq+(1:3*mr.nq))=0; e=1e-6;
            a=premelting_switch_ring_residual(z+e*direction,z,100,mr,bc,pr,f,zeros(mr.nq,1));
            b=premelting_switch_ring_residual(z-e*direction,z,100,mr,bc,pr,f,zeros(mr.nq,1));
            tc.verifyLessThan(norm((a-b)/(2*e)-J*direction)/norm(J*direction),1e-6);
        end
        function residualOnlyAndImmutableCachePreserveMixedPhaseEquations(tc)
            [~,m,~,p,f,z]=TestPremeltingResearch.fixture(); p.viscosityModel='iapws2008';
            n=m.nq; history=zeros(n,1); old=z;
            bc=thermo_interface_boundary(m,279.15,269.15,2e5,false);
            oldState=premelting_switch_ring_state(old,m,p,f,history);
            for temperature=[269.15,273.15,279.15]
                z(4*n+(1:n))=temperature;
                z(n+(1:3*n))=repmat([.3,.5,.2],n,1);
                [r,J]=premelting_switch_ring_residual(z,old,100,m,bc,p,f,history);
                value=premelting_switch_ring_residual(z,old,100,m,bc,p,f,history);
                [cached,Jcached]=premelting_switch_ring_residual(z,old,100,m,bc,p,f,history,[],oldState);
                tc.verifyEqual(value,r); tc.verifyEqual(cached,r); tc.verifyEqual(Jcached,J);
            end
        end
        function trivialPhaseEliminationPreservesNewtonEquations(tc)
            [~,m,~,p,f,z]=TestPremeltingResearch.fixture(); n=m.nq;
            bc=thermo_interface_boundary(m,279.15,269.15,2e5,false);
            for temperature=[269.15,273.15,279.15]
                z(4*n+(1:n))=temperature;
                [r,J]=premelting_switch_ring_residual(z,z,100,m,bc,p,f,zeros(n,1));
                old=premelting_newton_direction(J,-r,z(n+(1:3*n)),n,false);
                [fast,info]=premelting_newton_direction(J,-r,z(n+(1:3*n)),n,true);
                tc.verifyLessThan(norm(fast-old)/max(norm(old),1e-20),1e-8);
                tc.verifyLessThan(info.remainingRows,info.originalSize);
            end
        end
        function predictorUsesCommittedRateWithoutRelaxingStepLimit(tc)
            o=struct('stepController','transition_predictor','dtDays',1, ...
                'maximumPhaseIncrement',.01,'maximumPressureIncrement_Pa',2e4, ...
                'maximumDamageIncrement',.025,'maximumActiveIcePressureIncrement_Pa',1e4);
            tc.verifyEqual(premelting_step_proposal(100,.009,0,0,0,o),100,'AbsTol',1e-12);
            tc.verifyEqual(premelting_step_proposal(100,0,2e4,0,0,o),90,'AbsTol',1e-12);
            tc.verifyEqual(premelting_step_proposal(100,0,0,0,1e4,o),90,'AbsTol',1e-12);
            tc.verifyEqual(premelting_step_proposal(100,0,0,0,0,o),150);
            tc.verifyEqual(premelting_step_proposal(86400,0,0,0,0,o),86400);
            o.minimumStepSeconds=1;
            tc.verifyEqual(premelting_step_proposal(.5,0,2e4,0,0,o),1);
        end
    end
    methods (Static,Access=private)
        function removeTestFolder(folder,root)
            target=char(java.io.File(folder).getCanonicalPath());
            base=char(java.io.File(root).getCanonicalPath());
            assert(startsWith(target,[base,filesep]),'PhaseSwitch:TestCleanup','Test cleanup must remain within its workspace.');
            if isfolder(target), rmdir(target,'s'); end
        end
        function [m,mr,p,pr,f,z]=fixture()
            p=joint_threephase_parameters(); f=premelting_switch_parameters(); c=joint_cohesive_parameters(); c.enabled=false;
            m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
            [mr,pr]=reduce_premelting_cyclic_symmetry(m,p); n=mr.nq; ns=numel(mr.capacity);
            z=[.1*ones(n,1);ones(n,1);zeros(2*n,1);279.15*ones(n,1);.1*ones(n,1);279.15*ones(ns,1)];
        end
    end
end
