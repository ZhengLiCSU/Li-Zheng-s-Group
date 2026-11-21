classdef TestTHMBulkConstitutive < matlab.unittest.TestCase
    properties (TestParameter)
        PorosityFeedback=struct('off',false,'on',true)
    end
    methods (Test)
        function propertiesOnlyAndCacheMatchFullLaw(testCase,PorosityFeedback)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(PorosityFeedback);
            p.enthalpyReference_Jkg=1.234e5;
            full=thm_bulk_constitutive(T,h,e,p);
            base=frozen_soil_constitutive(T,h,p,false);
            value=thm_bulk_constitutive(T,h,e,p,[],false);
            cached=thm_bulk_constitutive(T,h,e,p,base,false);
            cachedFull=thm_bulk_constitutive(T,h,e,p,frozen_soil_constitutive(T,h,p),false);
            names=fieldnames(full);
            derivativeNames={'slT','slH','siT','siH','dMassT','dMassH', ...
                'dEnergyT','dEnergyH','dkT','dkH','dKT','dKH', ...
                'dMassEps','dEnergyEps','dkEps','dKEps', ...
                'dReferencePoreVolumeEps','dCurrentPorosityEps'};
            for j=1:numel(names)
                name=names{j};
                if ismember(name,derivativeNames)
                    testCase.verifyEmpty(value.(name),name);
                    testCase.verifyEmpty(cached.(name),name);
                    testCase.verifyEmpty(cachedFull.(name),name);
                else
                    testCase.verifyEqual(value.(name),full.(name),name,'AbsTol',0);
                    testCase.verifyEqual(cached.(name),full.(name),name,'AbsTol',0);
                    testCase.verifyEqual(cachedFull.(name),full.(name),name,'AbsTol',0);
                end
            end
            testCase.verifyError(@()thm_bulk_constitutive(T,h,e,p,base,true), ...
                'THMBulk:FrozenState');
        end

        function zeroStrainPreservesFrozenKernelAndPhasePartials(testCase)
            [T,h,~,p]=TestTHMBulkConstitutive.fixture(false);
            base=frozen_soil_constitutive(T,h,p);
            a=thm_bulk_constitutive(T,h,0,p);
            testCase.verifyEqual(a.mass,base.mass,'AbsTol',0);
            testCase.verifyEqual(a.energy,base.energy,'AbsTol',0);
            testCase.verifyEqual(a.K,base.K,'AbsTol',0);
            testCase.verifyEqual(a.k,base.k,'AbsTol',2e-15);
            testCase.verifyEqual(a.slT,base.slT,'AbsTol',0);
            testCase.verifyEqual(a.slH,base.slH,'AbsTol',0);
            testCase.verifyEqual(a.siT,base.siT,'AbsTol',0);
            testCase.verifyEqual(a.siH,base.siH,'AbsTol',0);
            testCase.verifyEqual(a.currentPorosity,p.porosity*ones(size(T)),'AbsTol',0);
        end

        function valuesMatchPriorAdditiveReferenceVolumeLaw(testCase)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(false);
            base=frozen_soil_constitutive(T,h,p); a=thm_bulk_constitutive(T,h,e,p);
            Vp=p.porosity+p.mechanics.biot*e; phi=Vp./(1+e);
            expectedMass=base.mass+(Vp-p.porosity).*(base.sl+base.si);
            expectedEnergy=base.energy+p.rho*(Vp-p.porosity).* ...
                (base.sl.*(p.L+p.cw*T)+base.si.*(p.ci*T));
            expectedK=(1-phi)*p.ks+phi.*(base.sl*p.kw+base.si*p.ki+base.sg*p.ka);
            testCase.verifyEqual(a.mass,expectedMass,'AbsTol',0);
            testCase.verifyEqual(a.energy,expectedEnergy,'AbsTol',0);
            testCase.verifyEqual(a.k,expectedK,'AbsTol',2e-15);
            testCase.verifyEqual(a.currentThetaW+a.currentThetaI+a.currentThetaA,phi,'AbsTol',2e-15);
        end

        function temperaturePartialsMatchFiniteDifferences(testCase,PorosityFeedback)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(PorosityFeedback); delta=2e-6;
            a=thm_bulk_constitutive(T,h,e,p);
            plus=thm_bulk_constitutive(T+delta,h,e,p); minus=thm_bulk_constitutive(T-delta,h,e,p);
            testCase.verifyEqual(a.dMassT,(plus.mass-minus.mass)/(2*delta),'AbsTol',2e-10,'RelTol',3e-6);
            testCase.verifyEqual(a.dEnergyT,(plus.energy-minus.energy)/(2*delta),'AbsTol',.1,'RelTol',3e-6);
            testCase.verifyEqual(a.dkT,(plus.k-minus.k)/(2*delta),'AbsTol',2e-9,'RelTol',3e-6);
            testCase.verifyEqual(a.dKT,(plus.K-minus.K)/(2*delta),'AbsTol',1e-13,'RelTol',3e-6);
        end

        function headPartialsMatchFiniteDifferences(testCase,PorosityFeedback)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(PorosityFeedback); delta=2e-6;
            a=thm_bulk_constitutive(T,h,e,p);
            plus=thm_bulk_constitutive(T,h+delta,e,p); minus=thm_bulk_constitutive(T,h-delta,e,p);
            testCase.verifyEqual(a.dMassH,(plus.mass-minus.mass)/(2*delta),'AbsTol',2e-10,'RelTol',3e-6);
            testCase.verifyEqual(a.dEnergyH,(plus.energy-minus.energy)/(2*delta),'AbsTol',.1,'RelTol',3e-6);
            testCase.verifyEqual(a.dkH,(plus.k-minus.k)/(2*delta),'AbsTol',2e-9,'RelTol',3e-6);
            testCase.verifyEqual(a.dKH,(plus.K-minus.K)/(2*delta),'AbsTol',1e-13,'RelTol',3e-6);
        end

        function strainPartialsMatchFiniteDifferences(testCase,PorosityFeedback)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(PorosityFeedback); delta=1e-7;
            a=thm_bulk_constitutive(T,h,e,p);
            plus=thm_bulk_constitutive(T,h,e+delta,p); minus=thm_bulk_constitutive(T,h,e-delta,p);
            testCase.verifyEqual(a.dMassEps,(plus.mass-minus.mass)/(2*delta),'AbsTol',2e-8,'RelTol',3e-6);
            testCase.verifyEqual(a.dEnergyEps,(plus.energy-minus.energy)/(2*delta),'AbsTol',.3,'RelTol',3e-6);
            testCase.verifyEqual(a.dkEps,(plus.k-minus.k)/(2*delta),'AbsTol',2e-8,'RelTol',3e-6);
            testCase.verifyEqual(a.dKEps,(plus.K-minus.K)/(2*delta),'AbsTol',1e-13,'RelTol',3e-6);
            testCase.verifyEqual(a.dCurrentPorosityEps, ...
                (plus.currentPorosity-minus.currentPorosity)/(2*delta),'AbsTol',2e-8,'RelTol',3e-6);
        end

        function enthalpyDatumRemainsConjugateToEveryMassPartial(testCase)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(true);
            a=thm_bulk_constitutive(T,h,e,p); datum=1.234e5; p.enthalpyReference_Jkg=datum;
            b=thm_bulk_constitutive(T,h,e,p);
            testCase.verifyEqual(b.energy-a.energy,p.rho*datum*a.mass,'AbsTol',5e-7,'RelTol',5e-14);
            testCase.verifyEqual(b.dEnergyT-a.dEnergyT,p.rho*datum*a.dMassT,'AbsTol',5e-7,'RelTol',5e-13);
            testCase.verifyEqual(b.dEnergyH-a.dEnergyH,p.rho*datum*a.dMassH,'AbsTol',5e-7,'RelTol',5e-13);
            testCase.verifyEqual(b.dEnergyEps-a.dEnergyEps,p.rho*datum*a.dMassEps,'AbsTol',5e-7,'RelTol',5e-13);
        end

        function optionalExistingFrozenStateMatchesFreshEvaluation(testCase)
            [T,h,e,p]=TestTHMBulkConstitutive.fixture(true);
            base=frozen_soil_constitutive(T,h,p);
            fresh=thm_bulk_constitutive(T,h,e,p);
            cached=thm_bulk_constitutive(T,h,e,p,base);
            testCase.verifyEqual(cached,fresh);
        end

        function kinematicVolumeAndBiotBoundsAreEnforced(testCase)
            [T,h,~,p]=TestTHMBulkConstitutive.fixture(false);
            testCase.verifyError(@()thm_bulk_constitutive(T,h,-1,p),'THMBulk:PoreVolumeDomain');
            p.mechanics.biot=1.01;
            testCase.verifyError(@()thm_bulk_constitutive(T,h,0,p),'THMBulk:Parameters');
        end
    end
    methods (Static,Access=private)
        function [T,h,e,p]=fixture(feedback)
            T=[5;3;-.8;-2]; h=[4;-2;3;-2]; e=[-.001;.003;.002;-.004];
            p=struct('rho',1000,'g',9.81,'L',3.34e5,'Tm',273.15, ...
                'porosity',.35,'thetaR',.05,'vgAlpha',.65,'vgN',1.67, ...
                'Ksat',1e-7,'Ss',1e-5,'iceImpedance',0, ...
                'cw',4180,'ci',2100,'rhoS',2650,'cs',800, ...
                'kw',.57,'ki',2.2,'ks',2.3,'ka',.025, ...
                'soilViscosity',struct('enabled',true,'referenceTemperature_C',20), ...
                'mechanics',struct('biot',.65), ...
                'thm',struct('permeabilityPorosityFeedback',feedback));
        end
    end
end
