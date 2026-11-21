classdef TestFrozenSoilConstitutive < matlab.unittest.TestCase
    % Independent checks for the pointwise equal-density TH law and tangents.
    methods (Test)
        function propertiesOnlyMatchesFullLawExactly(testCase)
            p=TestFrozenSoilConstitutive.fixture(); p.iceImpedance=1.3;
            p.enthalpyReference_Jkg=-1.234e5;
            T=[15 -1 -.025;-5 -.5 0]; h=[-2 -2 -100;3 -.25 0];
            for enabled=[false true]
                p.soilViscosity.enabled=enabled;
                full=frozen_soil_constitutive(T,h,p);
                value=frozen_soil_constitutive(T,h,p,false);
                names=fieldnames(full);
                for j=1:numel(names)
                    name=names{j};
                    if isempty(value.(name))
                        testCase.verifyTrue(ismember(name,{'slT','slH','siT','siH', ...
                            'dMassT','dMassH','dEnergyT','dEnergyH','dkT','dkH','dKT','dKH'}));
                    else
                        testCase.verifyEqual(value.(name),full.(name),name,'AbsTol',0);
                    end
                end
                testCase.verifyEmpty(value.dMassT);
                testCase.verifyEmpty(value.dKT);
            end
        end

        function singleOutputViscosityRetainsIdenticalValues(testCase)
            T=[-20 -.1 0 20 110];
            muOnly=water_viscosity_iapws(T);
            [mu,dmu]=water_viscosity_iapws(T);
            testCase.verifyEqual(muOnly,mu,'AbsTol',0);
            testCase.verifyTrue(all(isfinite(dmu)));
            testCase.verifyLessThan(dmu,zeros(size(T)));
        end

        function matchesUnweightedV14Formula(testCase)
            p=TestFrozenSoilConstitutive.fixture();
            T=[15;-1;-.025;-5;-.5;3]; h=[-2;-2;-100;3;-.25;5];
            actual=frozen_soil_constitutive(T,h,p);
            expected=TestFrozenSoilConstitutive.v14Formula(T,h,p);
            names=fieldnames(expected);
            for k=1:numel(names)
                name=names{k};
                testCase.verifyEqual(actual.(name),expected.(name), ...
                    ['v14 field: ' name],'AbsTol',1e-12,'RelTol',3e-13);
            end
        end

        function analyticTangentsMatchIndependentFiniteDifferences(testCase)
            p=TestFrozenSoilConstitutive.fixture(); p.iceImpedance=1.3;
            T=[5;-.04;-1;-5;-.1]; h=[-2;-100;-2;3;-.25];
            a=frozen_soil_constitutive(T,h,p); epsT=2e-7; epsH=2e-6;
            tp=frozen_soil_constitutive(T+epsT,h,p);
            tm=frozen_soil_constitutive(T-epsT,h,p);
            hp=frozen_soil_constitutive(T,h+epsH,p);
            hm=frozen_soil_constitutive(T,h-epsH,p);
            fields={'mass','energy','k','K'};
            derivT={'dMassT','dEnergyT','dkT','dKT'};
            derivH={'dMassH','dEnergyH','dkH','dKH'};
            absolute=[2e-9 2e-2 2e-8 2e-16];
            for k=1:numel(fields)
                f=fields{k};
                testCase.verifyEqual(a.(derivT{k}),(tp.(f)-tm.(f))/(2*epsT), ...
                    ['T tangent: ' f],'AbsTol',absolute(k),'RelTol',2e-5);
                testCase.verifyEqual(a.(derivH{k}),(hp.(f)-hm.(f))/(2*epsH), ...
                    ['h tangent: ' f],'AbsTol',absolute(k),'RelTol',2e-5);
            end
        end

        function phasePartitionAndApparentLiquidAreAdmissible(testCase)
            p=TestFrozenSoilConstitutive.fixture();
            [T,h]=ndgrid(linspace(-20,20,49),[-1e100 -1e5 -100 -2 -.01 0 5 1e100]);
            a=frozen_soil_constitutive(T,h,p);
            testCase.verifyGreaterThanOrEqual(a.sl,0);
            testCase.verifyGreaterThanOrEqual(a.si,0);
            testCase.verifyGreaterThanOrEqual(a.sg,-4e-16);
            testCase.verifyEqual(a.sl+a.si+a.sg,ones(size(T)),'AbsTol',3e-16);
            m=1-1/p.vgN; sr=p.thetaR/p.porosity;
            r=sr+(1-sr)*(1+(p.vgAlpha*max(-h,0)).^p.vgN).^(-m);
            testCase.verifyEqual(a.sl./(1-a.si),r,'AbsTol',2e-15);
            names=fieldnames(a);
            for k=1:numel(names)
                value=a.(names{k});
                testCase.verifyTrue(all(isfinite(value(:))),names{k});
            end
        end

        function enthalpyReferenceShiftIsConjugateToTotalWater(testCase)
            p=TestFrozenSoilConstitutive.fixture();
            T=[5;-.04;-1;-5;-.1]; h=[-2;-100;-2;3;-.25];
            a=frozen_soil_constitutive(T,h,p);
            p.enthalpyReference_Jkg=-2.5e5;
            b=frozen_soil_constitutive(T,h,p);
            scale=p.rho*p.enthalpyReference_Jkg;
            testCase.verifyEqual(b.energy-a.energy,scale*a.mass, ...
                'AbsTol',5e-8,'RelTol',2e-14);
            testCase.verifyEqual(b.dEnergyT-a.dEnergyT,scale*a.dMassT, ...
                'AbsTol',5e-8,'RelTol',2e-14);
            testCase.verifyEqual(b.dEnergyH-a.dEnergyH,scale*a.dMassH, ...
                'AbsTol',5e-8,'RelTol',2e-14);
            testCase.verifyEqual(b.mass,a.mass,'AbsTol',0);
            testCase.verifyEqual(b.K,a.K,'AbsTol',0);
            p.enthalpyReference_Jkg=0;
            testCase.verifyEqual(frozen_soil_constitutive(T,h,p),a);
        end

        function saturatedEqualDensityMassHasNoThermalSource(testCase)
            p=TestFrozenSoilConstitutive.fixture();
            T=[-10;-.5;-.01;1;12]; a=frozen_soil_constitutive(T,4,p);
            testCase.verifyEqual(a.mass,repmat(p.porosity+4*p.Ss,size(T)), ...
                'AbsTol',2e-16);
            testCase.verifyEqual(a.dMassT,zeros(size(T)),'AbsTol',2e-15);
            testCase.verifyEqual(a.dMassH,repmat(p.Ss,size(T)),'AbsTol',0);
        end

        function warmRetentionDoesNotContainThermalHead(testCase)
            p=TestFrozenSoilConstitutive.fixture(); h=[-100;-2;.1;5];
            a=frozen_soil_constitutive(2,h,p); b=frozen_soil_constitutive(12,h,p);
            testCase.verifyEqual(a.sl,b.sl,'AbsTol',0);
            testCase.verifyEqual(a.mass,b.mass,'AbsTol',0);
            testCase.verifyEqual(a.si,zeros(size(h)),'AbsTol',0);
            testCase.verifyEqual(a.dMassT,zeros(size(h)),'AbsTol',0);
        end

        function viscosityTangentAndDisabledMode(testCase)
            p=TestFrozenSoilConstitutive.fixture(); T=[-2;5;20];
            a=frozen_soil_constitutive(T,3,p);
            [mu,dmu]=water_viscosity_iapws(T);
            % At fixed positive h, the warm derivative is pure viscosity.
            testCase.verifyEqual(a.dKT(2:3),-a.K(2:3).*dmu(2:3)./mu(2:3), ...
                'AbsTol',1e-22,'RelTol',2e-14);
            p.soilViscosity.enabled=false;
            b=frozen_soil_constitutive(T,3,p);
            testCase.verifyEqual(b.dKT(2:3),zeros(2,1),'AbsTol',0);
            testCase.verifyEqual(b.K(2:3),repmat(p.Ksat,2,1),'AbsTol',1e-22);
        end

        function flatBranchTangentsAtSwitches(testCase)
            p=TestFrozenSoilConstitutive.fixture();
            a=frozen_soil_constitutive(0,0,p);
            testCase.verifyEqual(a.dMassT,0,'AbsTol',0);
            testCase.verifyEqual(a.dMassH,0,'AbsTol',0);
            testCase.verifyEqual(a.dkT,0,'AbsTol',0);
            testCase.verifyEqual(a.dkH,0,'AbsTol',0);
            testCase.verifyTrue(isfinite(a.dKT));
            testCase.verifyEqual(a.dKH,0,'AbsTol',0);
            % Equality F=r chooses the unfrozen min branch. T<0 avoids T=0.
            T=-.1; h=p.L/(p.g*p.Tm)*T;
            b=frozen_soil_constitutive(T,h,p);
            testCase.verifyEqual(b.si,0,'AbsTol',0);
            testCase.verifyEqual(b.dMassT,0,'AbsTol',0);
            testCase.verifyGreaterThan(b.dMassH,0);
        end

        function zeroPermeabilityHasZeroTangents(testCase)
            p=TestFrozenSoilConstitutive.fixture(); p.Ksat=0;
            a=frozen_soil_constitutive([-2 5],[-1 3],p);
            testCase.verifyEqual(a.K,[0 0],'AbsTol',0);
            testCase.verifyEqual(a.dKT,[0 0],'AbsTol',0);
            testCase.verifyEqual(a.dKH,[0 0],'AbsTol',0);
        end

        function scalarExpansionAndInvalidInputs(testCase)
            p=TestFrozenSoilConstitutive.fixture();
            a=frozen_soil_constitutive(-2,[-3 5],p);
            b=frozen_soil_constitutive([-2 -2],[-3 5],p);
            testCase.verifyEqual(a,b);
            testCase.verifyError(@()frozen_soil_constitutive([1 2],[1;2],p), ...
                'frozen_soil_constitutive:SizeMismatch');
            testCase.verifyError(@()frozen_soil_constitutive(NaN,0,p), ...
                'frozen_soil_constitutive:InvalidState');
            p.thetaR=p.porosity;
            testCase.verifyError(@()frozen_soil_constitutive(1,0,p), ...
                'frozen_soil_constitutive:InvalidParameters');
        end
    end

    methods (Static,Access=private)
        function p=fixture()
            p=struct('rho',1000,'g',9.81,'L',3.34e5,'Tm',273.15, ...
                'porosity',.35,'thetaR',.05,'vgAlpha',.65,'vgN',1.67, ...
                'Ksat',9.81e-8,'Ss',1e-5,'iceImpedance',0, ...
                'cw',4180,'ci',2100,'rhoS',2650,'cs',800, ...
                'kw',.57,'ki',2.2,'ks',2.3,'ka',.025, ...
                'soilViscosity',struct('enabled',true,'referenceTemperature_C',20));
        end

        function a=v14Formula(T,h,p)
            % Independent direct transcription of v14's unweighted state.
            m=1-1/p.vgN; sr=p.thetaR/p.porosity;
            R=@(x)sr+(1-sr)*(1+(p.vgAlpha*max(-x,0)).^p.vgN).^(-m);
            r=R(h); F=R(-p.L/(p.g*p.Tm)*max(-T,0));
            a.sl=min(r,F); a.si=max(0,1-F./r); a.sg=1-a.sl-a.si;
            a.thetaW=p.porosity*a.sl; a.thetaI=p.porosity*a.si;
            a.thetaA=p.porosity*a.sg; a.total=a.thetaW+a.thetaI;
            a.elasticWater=p.Ss*max(h,0); a.mass=a.total+a.elasticWater;
            a.energy=(1-p.porosity)*p.rhoS*p.cs*T+p.rho*( ...
                a.thetaW.*(p.L+p.cw*T)+a.thetaI.*p.ci.*T+ ...
                a.elasticWater.*(p.L+p.cw*T));
            a.k=(1-p.porosity)*p.ks+a.thetaW*p.kw+a.thetaI*p.ki+a.thetaA*p.ka;
            Se=max(0,min(1,(a.sl-sr)/(1-sr)));
            kr=sqrt(Se).*(1-(1-Se.^(1/m)).^m).^2;
            a.K=p.Ksat*kr.*10.^(-p.iceImpedance*a.si).* ...
                (water_viscosity_iapws(p.soilViscosity.referenceTemperature_C)./ ...
                water_viscosity_iapws(T));
        end
    end
end
