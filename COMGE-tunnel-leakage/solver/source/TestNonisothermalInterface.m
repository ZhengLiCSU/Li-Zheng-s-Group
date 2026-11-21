classdef TestNonisothermalInterface < matlab.unittest.TestCase
    % Independent face identities, analytic limits, coupled conservation.
    methods (Test)
        function temperatureGradientDoesNotCreateMassFlux(testCase)
            [f,~,a]=joint_nonisothermal_face([2e5,270],[2e5,290],0.01,2);
            testCase.verifyEqual(f,[0;-40],'AbsTol',1e-14);
            testCase.verifyGreaterThan(a.entropyProduction_WmK,0);
        end
        function faceEntropyEqualsThermodynamicIdentity(testCase)
            [~,~,a]=joint_nonisothermal_face([3e6,269],[2e5,286],0.02,0.3);
            testCase.verifyEqual(a.entropyProduction_WmK,a.directEntropyProduction_WmK,'RelTol',1e-12);
            testCase.verifyGreaterThanOrEqual([a.hydraulic_WmK,a.fourier_WmK,a.upwind_WmK],zeros(1,3));
        end
        function reverseFlowProducesEntropy(testCase)
            [~,~,a]=joint_nonisothermal_face([0,290],[3e6,269],0.02,0.3);
            testCase.verifyLessThan(a.massFlux_kgms,0);
            testCase.verifyEqual(a.entropyProduction_WmK,a.directEntropyProduction_WmK,'RelTol',1e-12);
        end
        function faceTangentMatchesIndependentDifference(testCase)
            x=[3e5,276,1e5,274,1e-4,0.01]; scale=[1e6,1,1e6,1,1e-4,0.01];
            [~,J]=joint_nonisothermal_face(x(1:2),x(3:4),x(5),x(6));
            fd=TestNonisothermalInterface.faceDifference(x,scale);
            testCase.verifyEqual(J.*scale,fd,'RelTol',2e-6,'AbsTol',1e-7);
        end
        function allIceBlocksMassButConductsHeat(testCase)
            [f,~,a]=joint_nonisothermal_face([2e5,268],[0,270],0,2);
            testCase.verifyEqual(f,[0;-4]);
            testCase.verifyGreaterThan(a.fourier_WmK,0);
        end
        function ringUsesIndependentFacesAndReciprocalPorts(testCase)
            [m,p,z]=TestNonisothermalInterface.fixture(); s=thermo_interface_state(z,m,p);
            virtual=sin((1:numel(s.u))')*1e-7;
            load=m.gapOperator'*(m.weight*p.wetAreaFraction.*s.pressure);
            testCase.verifyEqual(load'*virtual, ...
                sum(m.weight*p.wetAreaFraction.*s.pressure.*(m.gapOperator*virtual)),'AbsTol',1e-15);
            testCase.verifyEqual(numel(unique(m.mesh.connectedComponentID)),8);
        end
        function coupledTangentMatchesDirectionalDifference(testCase)
            [m,p,z]=TestNonisothermalInterface.fixture(); n=numel(m.weight);
            z(n+(1:n))=0.2; z(2*n+(1:n))=273.14;
            bc=thermo_interface_boundary(m,276,274);
            [~,J]=thermo_interface_residual(z,z,3600,m,bc,p);
            direction=[1e5*cos((1:n)');0.05*sin((1:n)');cos((1:numel(z)-2*n)')];
            h=1e-5; fd=(thermo_interface_residual(z+h*direction,z,3600,m,bc,p)- ...
                thermo_interface_residual(z-h*direction,z,3600,m,bc,p))/(2*h);
            testCase.verifyEqual(J*direction,fd,'RelTol',2e-5,'AbsTol',2e-6);
        end
        function warmFlowBalancesLocalMassAndTotalEnergy(testCase)
            [m,p,z]=TestNonisothermalInterface.fixture(); bc=thermo_interface_boundary(m,279.15,279.15);
            [~,a]=thermo_interface_step(z,60,m,bc,p);
            testCase.verifyLessThan(max(abs(a.localMassDefect_kgm)),1e-12);
            testCase.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
            testCase.verifyLessThan(abs(sum(a.jointHeat_Wm)+sum(a.solidHeat_Wm)),1e-10);
            testCase.verifyGreaterThan(a.entropyProduced_JmK,0);
        end
        function isolatedRingConservesEnergyDuringThermalRelaxation(testCase)
            [m,p,z]=TestNonisothermalInterface.fixture(); n=numel(m.weight);
            z(2*n+(1:n))=276; m.bathConductance(:)=0;
            bc=thermo_interface_boundary(m,279.15,279.15,0,true); bc.insulatedHeat(:)=true;
            [~,a]=thermo_interface_step(z,60,m,bc,p);
            testCase.verifyEqual(a.boundaryEnergy_Jm,0);
            testCase.verifyLessThan(abs(a.totalEnergyDefect_Jm),a.energyAcceptanceBound_Jm);
            testCase.verifyGreaterThan(a.entropyProduced_JmK,0);
        end
        function noFlowThermalEquilibriumRemainsStationary(testCase)
            [m,p,z]=TestNonisothermalInterface.fixture();
            bc=thermo_interface_boundary(m,279.15,279.15,0,true);
            [zn,a]=thermo_interface_step(z,1000,m,bc,p);
            testCase.verifyEqual(zn,z,'AbsTol',1e-9);
            testCase.verifyLessThan(abs(a.totalEnergyDefect_Jm),1e-6);
        end
        function latentHeatMatchesIndependentAdiabaticRoot(testCase)
            p=joint_thermodynamic_parameters(); m=TestNonisothermalInterface.singleCell(p);
            old=[0;0;272.15;272.15]; bc=struct('reservoirPressure_Pa',[0;0], ...
                'reservoirTemperature_K',[272.15;272.15],'sealedMass',[true;true], ...
                'insulatedHeat',[true;true],'outletSeepage',[false;false],'solidBathTemperature_K',272.15);
            [z,a]=thermo_interface_step(old,1,m,bc,p);
            [expected,~]=joint_adiabatic_equilibrium([0,0.005,0,272.15],p);
            testCase.verifyEqual(z(3),expected(4),'AbsTol',1e-8);
            testCase.verifyEqual(a.state.x(3),expected(3),'AbsTol',1e-10);
            testCase.verifyLessThan(abs(a.totalEnergyDefect_Jm),1e-8);
        end
        function frozenTensileLimitIsReportedWithoutClipping(testCase)
            [m,p,z]=TestNonisothermalInterface.fixture(); n=numel(m.weight);
            z(1:n)=p.minimumPressure_Pa+1; z(n+(1:n))=1;
            a=diagnose_thermo_interface_domain(z,m,p);
            testCase.verifyEqual(a.classification,'saturated_single_pressure_domain_exhausted');
            testCase.verifyEqual(numel(a.fullyFrozenPressureLimitedCells),n);
            testCase.verifyEqual(a.minimumPressure_Pa,p.minimumPressure_Pa+1);
        end
        function publicRunnerDefaultsResolvePressureChanges(testCase)
            [m,p,~]=TestNonisothermalInterface.fixture(); n=numel(m.weight);
            out=solve_nonisothermal_ring(m,p,struct('cycles',0.01,'temperatureAmplitude_K',0,'verbose',false));
            testCase.verifyTrue(out.completed);
            testCase.verifyEqual(out.options.maximumPhaseIncrement,0.025);
            testCase.verifyLessThanOrEqual(max(abs(diff(out.states(1:n,:),1,2)),[],'all'),1e4);
            testCase.verifyEqual(out.time_s(end),0.01*20*86400,'AbsTol',1e-6);
        end
    end
    methods (Static, Access=private)
        function [m,p,z]=fixture()
            p=joint_thermodynamic_parameters();
            m=build_thermo_interface_ring(make_thm_verification_mesh(),2,p);
            n=numel(m.weight); ns=numel(m.capacity);
            z=[1e5*ones(n,1);zeros(n,1);279.15*ones(n+ns,1)];
        end
        function fd=faceDifference(x,scale)
            fd=zeros(2,6); h=1e-5;
            for j=1:6
                a=x; b=x; a(j)=a(j)+h*scale(j); b(j)=b(j)-h*scale(j);
                fa=joint_nonisothermal_face(a(1:2),a(3:4),a(5),a(6));
                fb=joint_nonisothermal_face(b(1:2),b(3:4),b(5),b(6));
                fd(:,j)=(fa-fb)/(2*h);
            end
        end
        function m=singleCell(p)
            m=struct('weight',1,'capacity',1,'stiffness',p.normalStiffness_Pam, ...
                'thermalStiffness',0,'thermalLoad',0,'gapOperator',1, ...
                'displacementPressure',p.wetAreaFraction/p.normalStiffness_Pam, ...
                'displacementTemperature',0,'heatExchange',0,'exchangeRowSum',0, ...
                'exchangeColumnSum',0,'bathConductance',0,'referenceTemperature_K',272.15, ...
                'conductivityWater',0.6,'conductivityIce',2.2,'faces',[0,1;1,0], ...
                'faceDistance',[0.5;0.5],'boundaryID',[1;2],'divergence',[-1,1]);
        end
    end
end
