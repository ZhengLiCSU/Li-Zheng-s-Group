classdef TestJointThermodynamics < matlab.unittest.TestCase
    % Analytical identities, physical limiting cases and conservative coupling.
    methods (Test)
        function phaseVolumesFillExactlyAvailableAperture(testCase)
            [x,p] = TestJointThermodynamics.mixedState();
            s = joint_helmholtz_state(x,p);
            testCase.verifyEqual(s.liquidVolume_m+s.iceVolume_m,s.aperture_m,'AbsTol',1e-20);
        end
        function pressureNormalWorkUsesSameVolumeDerivative(testCase)
            [x,p] = TestJointThermodynamics.mixedState();
            p.wetAreaFraction = 0.25;
            s = joint_helmholtz_state(x,p,false);
            testCase.verifyEqual(s.fluidNormalTraction_Pa,-0.25*s.pressure_Pa);
        end
        function helmholtzGradientMatchesFiniteDifferences(testCase)
            [x,p] = TestJointThermodynamics.mixedState();
            s = joint_helmholtz_state(x,p);
            numerical = TestJointThermodynamics.energyGradient(x,p);
            scale = [1e-6,0.005,0.005,1];
            testCase.verifyEqual(s.gradient.*scale,numerical,'AbsTol',2e-7);
        end
        function exactHessianMatchesFiniteDifferences(testCase)
            [x,p] = TestJointThermodynamics.mixedState();
            s = joint_helmholtz_state(x,p);
            numerical = TestJointThermodynamics.energyHessian(x,p);
            scale = diag([1e-6,0.005,0.005,1]);
            testCase.verifyEqual(scale*s.hessian*scale,numerical,'RelTol',2e-8,'AbsTol',2e-8);
        end
        function maxwellSymmetryAndIsothermalConvexity(testCase)
            [x,p] = TestJointThermodynamics.mixedState();
            s = joint_helmholtz_state(x,p);
            scale = diag([1e-6,0.005,0.005,1]); h = scale*s.hessian*scale;
            testCase.verifyEqual(h,h','RelTol',1e-14,'AbsTol',1e-12);
            testCase.verifyGreaterThanOrEqual(min(eig(h(1:3,1:3))),-1e-12);
        end
        function temperatureCurvatureGivesPositiveHeatCapacity(testCase)
            [x,p] = TestJointThermodynamics.mixedState(); s = joint_helmholtz_state(x,p);
            testCase.verifyEqual(-x(4)*s.hessian(4,4),s.heatCapacity_Jm2K,'RelTol',1e-14);
        end
        function internalEnergyUsesLegendreTransform(testCase)
            [x,p] = TestJointThermodynamics.mixedState(); p.externalTraction_Pa = 1e5;
            s = joint_helmholtz_state(x,p);
            testCase.verifyEqual(s.internalEnergy_Jm2, ...
                s.freeEnergy_Jm2+x(4)*s.entropy_Jm2K+p.externalTraction_Pa*x(1),'AbsTol',2e-12);
        end
        function pressureMeltingSlopeMatchesClapeyron(testCase)
            p = joint_thermodynamic_parameters();
            slope = TestJointThermodynamics.meltingTemperature(100,p)-p.meltingTemperature_K;
            analytical = -p.meltingTemperature_K*(1/p.rhoIce-1/p.rhoWater)/p.latentHeat_Jkg;
            testCase.verifyEqual(slope/100,analytical,'RelTol',2e-6);
        end
        function bulkCoexistenceAtZeroPressure(testCase)
            p = joint_thermodynamic_parameters(); m = p.rhoWater*p.aperture0_m; i = 0.4*m;
            g = ((m-i)/p.rhoWater+i/p.rhoIce-p.aperture0_m)/p.wetAreaFraction;
            s = joint_helmholtz_state([g,m,i,p.meltingTemperature_K],p);
            testCase.verifyEqual(s.phaseAffinity_Jkg,0,'AbsTol',1e-9);
        end
        function iceFreeExpansionIsDensityRatio(testCase)
            p = joint_thermodynamic_parameters(); m = p.rhoWater*p.aperture0_m;
            g = m/p.rhoIce-p.aperture0_m;
            s = joint_helmholtz_state([g,m,m,272.15],p);
            testCase.verifyEqual(g/p.aperture0_m,p.rhoWater/p.rhoIce-1,'AbsTol',1e-15);
            testCase.verifyEqual(s.pressure_Pa,0,'AbsTol',1e-5);
        end
        function sealedFreezingConservesMassAndOpensWalls(testCase)
            p = joint_thermodynamic_parameters(); x = [0,p.rhoWater*p.aperture0_m,0,274.15];
            [next,info] = joint_thermodynamic_step(x,272.15,3600, ...
                struct('sealedEnds',[true;true],'tolerance',1e-12),p);
            testCase.verifyEqual(next(2),x(2),'AbsTol',1e-14);
            testCase.verifyEqual(next(3),next(2),'AbsTol',1e-14);
            testCase.verifyGreaterThan(next(1),0);
            testCase.verifyEqual(p.normalStiffness_Pam*next(1),info.state.pressure_Pa,'RelTol',1e-9);
            testCase.verifyGreaterThan(info.entropyProduced_JmK,0);
        end
        function rigidFreezingProducesClapeyronPressure(testCase)
            p = joint_thermodynamic_parameters(struct('normalStiffness_Pam',1e20));
            x = [0,p.rhoWater*p.aperture0_m,0,274.15];
            [next,info] = joint_thermodynamic_step(x,272.15,3600,struct('sealedEnds',[true;true]),p);
            reference = fzero(@(q)TestJointThermodynamics.affinityAtPressure(q,272.15,p),[0,5e7]);
            testCase.verifyEqual(info.state.pressure_Pa,reference,'RelTol',1e-7);
            testCase.verifyGreaterThan(next(3),0);
            testCase.verifyLessThan(next(3),next(2));
        end
        function closedFreezeThawReturnsInitialVolume(testCase)
            p = joint_thermodynamic_parameters(); x = [0,p.rhoWater*p.aperture0_m,0,274.15];
            frozen = joint_thermodynamic_step(x,272.15,3600,struct('sealedEnds',[true;true]),p);
            thawed = joint_thermodynamic_step(frozen,274.15,3600,struct('sealedEnds',[true;true]),p);
            testCase.verifyEqual(thawed,x,'AbsTol',1e-13);
        end
        function saturatedBranchRejectsCavitation(testCase)
            p = joint_thermodynamic_parameters(); x = [1e-6,p.rhoWater*p.aperture0_m,0,274.15];
            testCase.verifyError(@()joint_helmholtz_state(x,p),'JointThermo:CavitationDomain');
        end
        function smallCompressionDomainIsNotClipped(testCase)
            p = joint_thermodynamic_parameters(); x = [-4e-6,p.rhoWater*p.aperture0_m,0,274.15];
            testCase.verifyError(@()joint_helmholtz_state(x,p),'JointThermo:CompressionDomain');
        end
        function negativeIceMassIsRejected(testCase)
            testCase.verifyError(@()joint_helmholtz_state([0,0.005,-1e-5,274.15]),'JointThermo:StateDomain');
        end
        function unknownParametersAreRejected(testCase)
            testCase.verifyError(@()joint_thermodynamic_parameters(struct('fakeDamage',10)), ...
                'JointThermo:UnknownParameter');
        end
        function coupledResidualHasConsistentTangent(testCase)
            [x,p] = TestJointThermodynamics.mixedState(); old = repmat(x,2,1);
            xx = old(:,1:3); settings = TestJointThermodynamics.settings();
            [~,jacobian] = joint_interface_residual(xx,old,x(4),1,settings,p);
            [scaled,numerical] = TestJointThermodynamics.residualTangent(xx,old,x(4),settings,p,jacobian);
            testCase.verifyEqual(scaled,numerical,'RelTol',5e-7,'AbsTol',1e-7);
        end
        function hydraulicFluxIsZeroWithNoLiquidLayer(testCase)
            p = joint_thermodynamic_parameters(); m = p.rhoWater*p.aperture0_m;
            x = [m/p.rhoIce-p.aperture0_m,m,m,272.15];
            [~,~,flow] = joint_interface_residual(x(1:3),x,x(4),1,TestJointThermodynamics.settings(),p);
            testCase.verifyEqual(flow.massFlux_kgms,zeros(2,1));
        end
        function warmHydraulicsHaveConservativeMassAndPositiveDissipation(testCase)
            out = solve_thermodynamic_joint(struct('cells',4,'duration_s',3600,'dt_s',600, ...
                'initialTemperature_K',274.15,'temperatureFunction',@(t)274.15));
            testCase.verifyLessThan(abs(out.massLedger_kgm(4)),1e-10);
            testCase.verifyGreaterThan(out.outletReferenceVolume_Lm,0);
            testCase.verifyGreaterThanOrEqual(min(out.audits(:,5:8),[],'all'),-1e-8);
        end
        function primaryPressureMustSatisfyVolumeClosure(testCase)
            [x,p] = TestJointThermodynamics.mixedState();
            testCase.verifyError(@()joint_helmholtz_state(x,p,true,2e6), ...
                'JointThermo:PressureClosure');
        end
        function pressureCoordinateTangentMatchesFiniteDifferences(testCase)
            [x,p] = TestJointThermodynamics.mixedState(); old = repmat(x,2,1);
            z = repmat([x(1),1e6,x(3)/x(2)],2,1); settings = TestJointThermodynamics.settings();
            [~,jacobian] = joint_interface_primary_residual(z,old,x(4),1,settings,p);
            [scaled,numerical] = TestJointThermodynamics.primaryTangent(z,old,x(4),settings,p,jacobian);
            testCase.verifyEqual(scaled,numerical,'RelTol',5e-7,'AbsTol',1e-7);
        end
        function analyticalPoiseuilleLimitIsRecovered(testCase)
            out = solve_thermodynamic_joint(struct('parameters', ...
                struct('normalStiffness_Pam',1e20,'bulkWater_Pa',1e15), ...
                'cells',4,'duration_s',3600,'dt_s',600, ...
                'initialTemperature_K',274.15,'temperatureFunction',@(t)274.15));
            p = out.parameters;
            exact = p.aperture0_m^3/(12*p.viscosity_Pas*0.45)*2e5*3600*1000;
            testCase.verifyEqual(out.outletReferenceVolume_Lm,exact,'RelTol',1e-7);
        end
        function adiabaticFreezingReleasesLatentHeat(testCase)
            p = joint_thermodynamic_parameters(); x = [0,p.rhoWater*p.aperture0_m,0,272.15];
            [next,info] = joint_adiabatic_equilibrium(x,p);
            testCase.verifyGreaterThan(next(4),x(4));
            testCase.verifyLessThan(next(4),p.meltingTemperature_K);
            testCase.verifyGreaterThan(next(3),0);
            testCase.verifyLessThan(next(3)/next(2),0.02);
            testCase.verifyLessThan(abs(info.energyDefect_Jm2),1e-8);
            testCase.verifyGreaterThan(info.entropyProduction_Jm2K,0);
        end
        function reversibleGeometryDoesNotCreateDamage(testCase)
            p = joint_thermodynamic_parameters(); x = [0,p.rhoWater*p.aperture0_m,0,274.15];
            cold = joint_closed_equilibrium(x(2),272.15,p);
            warm = joint_closed_equilibrium(cold(2),274.15,p);
            testCase.verifyEqual(warm,x,'AbsTol',1e-13);
        end
        function finiteRatePhaseChangeHasNonnegativeDissipation(testCase)
            p = joint_thermodynamic_parameters(struct('phaseViscosity_Jm2skg2',1e10));
            x = [0,p.rhoWater*p.aperture0_m,0,274.15];
            [next,info] = joint_thermodynamic_step(x,272.15,3600,struct('sealedEnds',[true;true]),p);
            exact = 0.45*p.phaseViscosity_Jm2skg2*(next(3)-x(3))^2/3600;
            testCase.verifyEqual(info.phaseDissipation_Jm,exact,'RelTol',1e-8);
            testCase.verifyGreaterThan(info.entropyProduced_JmK,0);
            testCase.verifyLessThan(next(3),next(2));
        end
        function phaseSwitchRecoveryPreservesCompleteTwoYearLedger(testCase)
            out = solve_thermodynamic_joint(struct('cells',16,'duration_s',2*365*86400, ...
                'dt_s',21600,'initialTemperature_K',285.15));
            testCase.verifyEqual(out.time_s(end),2*365*86400);
            testCase.verifyLessThan(abs(out.massLedger_kgm(4)),1e-10);
            testCase.verifyLessThan(max(abs(out.audits(:,4))),1e-7);
            testCase.verifyEqual(max(out.stateHistory(:,3,:)./out.stateHistory(:,2,:),[],'all'),1,'AbsTol',1e-13);
            testCase.verifyGreaterThan(out.outletReferenceVolume_Lm,0);
        end
        function wetAreaScalesInventoryTractionAndFlowConsistently(testCase)
            out = solve_thermodynamic_joint(struct('parameters', ...
                struct('wetAreaFraction',0.25,'aperture0_m',1.25e-6, ...
                'normalStiffness_Pam',1e20,'bulkWater_Pa',1e15), ...
                'cells',4,'duration_s',3600,'dt_s',600, ...
                'initialTemperature_K',274.15,'temperatureFunction',@(t)274.15));
            exact = 0.25*(5e-6)^3/(12*out.parameters.viscosity_Pas*0.45)*2e5*3600*1000;
            testCase.verifyEqual(out.outletReferenceVolume_Lm,exact,'RelTol',1e-7);
            testCase.verifyLessThan(abs(out.massLedger_kgm(4)),1e-10);
        end
    end
    methods (Static, Access=private)
        function [x,p] = mixedState()
            p = joint_thermodynamic_parameters(); m = 0.005; i = 0.002; pressure = 1e6;
            v = (m-i)/p.rhoWater+i/p.rhoIce;
            c = (m-i)/(p.rhoWater*p.bulkWater_Pa)+i/(p.rhoIce*p.bulkIce_Pa);
            g = (v-c*pressure-p.aperture0_m)/p.wetAreaFraction;
            x = [g,m,i,273.1];
        end
        function grad = energyGradient(x,p)
            scale = [1e-6,0.005,0.005,1]; epsilon = 1e-6; grad = zeros(1,4);
            for k = 1:4
                d = zeros(1,4); d(k) = epsilon*scale(k);
                plus = joint_helmholtz_state(x+d,p); minus = joint_helmholtz_state(x-d,p);
                grad(k) = (plus.freeEnergy_Jm2-minus.freeEnergy_Jm2)/(2*epsilon);
            end
        end
        function h = energyHessian(x,p)
            scale = [1e-6,0.005,0.005,1]; epsilon = 1e-6; h = zeros(4);
            for k = 1:4
                d = zeros(1,4); d(k) = epsilon*scale(k);
                plus = joint_helmholtz_state(x+d,p); minus = joint_helmholtz_state(x-d,p);
                h(:,k) = ((plus.gradient-minus.gradient).*scale)'/(2*epsilon);
            end
        end
        function temperature = meltingTemperature(pressure,p)
            temperature = fzero(@(t)TestJointThermodynamics.affinityAtPressure(pressure,t,p), ...
                p.meltingTemperature_K);
        end
        function affinity = affinityAtPressure(pressure,temperature,p)
            thermal = temperature-p.meltingTemperature_K-temperature*log(temperature/p.meltingTemperature_K);
            affinity = (p.heatIce_JkgK-p.heatWater_JkgK)*thermal+ ...
                p.latentHeat_Jkg*(temperature/p.meltingTemperature_K-1)+p.capillaryShift_Jkg+ ...
                (1/p.rhoIce-1/p.rhoWater)*pressure-0.5*( ...
                1/(p.rhoIce*p.bulkIce_Pa)-1/(p.rhoWater*p.bulkWater_Pa))*pressure^2;
        end
        function settings = settings()
            settings = struct('length_m',0.45,'reservoirPressure_Pa',[2e5;0], ...
                'sealedEnds',[false;false],'outletSeepage',true);
        end
        function [scaled,numerical] = residualTangent(x,old,t,settings,p,jacobian)
            scale = repmat([1e-6;0.005;0.005],size(x,1),1); epsilon = 1e-6;
            numerical = zeros(numel(x),numel(x)); vector = reshape(x',[],1);
            for k = 1:numel(x)
                d = zeros(size(vector)); d(k) = epsilon*scale(k);
                plus = joint_interface_residual(reshape(vector+d,3,[])',old,t,1,settings,p);
                minus = joint_interface_residual(reshape(vector-d,3,[])',old,t,1,settings,p);
                numerical(:,k) = (plus-minus)/(2*epsilon);
            end
            scaled = full(jacobian*diag(scale));
        end
        function [scaled,numerical] = primaryTangent(z,old,t,settings,p,jacobian)
            scale = repmat([1e-6;1e6;1],size(z,1),1); epsilon = 1e-6;
            numerical = zeros(numel(z),numel(z)); vector = reshape(z',[],1);
            for k = 1:numel(z)
                d = zeros(size(vector)); d(k) = epsilon*scale(k);
                plus = joint_interface_primary_residual(reshape(vector+d,3,[])',old,t,1,settings,p);
                minus = joint_interface_primary_residual(reshape(vector-d,3,[])',old,t,1,settings,p);
                numerical(:,k) = (plus-minus)/(2*epsilon);
            end
            scaled = full(jacobian*diag(scale));
        end
    end
end
