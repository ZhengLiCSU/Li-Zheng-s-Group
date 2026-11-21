classdef TestFrozenSoilPropertiesJIT < matlab.unittest.TestCase
    % Public numeric kernel parity. Build once separately to enable MEX tests.
    properties (TestParameter)
        viscosityEnabled=struct('withViscosity',true,'withoutViscosity',false)
        referenceEnthalpy=struct('zeroDatum',0,'shiftedDatum',-1.234e5)
    end

    methods (Test)
        function numericEntryMatchesEveryProperty(testCase,viscosityEnabled,referenceEnthalpy)
            p=testCase.fixture(); p.soilViscosity.enabled=viscosityEnabled;
            p.enthalpyReference_Jkg=referenceEnthalpy;
            [T,h]=testCase.grid(); [v,schema]=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyProperties(actual,expected);
        end

        function compiledEntryMatchesEveryProperty(testCase,viscosityEnabled,referenceEnthalpy)
            testCase.assumeEqual(exist('frozen_soil_properties_jit_mex','file'),3, ...
                'Run build_frozen_soil_jit before this optional compiled check.');
            p=testCase.fixture(); p.soilViscosity.enabled=viscosityEnabled;
            p.enthalpyReference_Jkg=referenceEnthalpy;
            [T,h]=testCase.grid(); [v,schema]=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit_mex(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyProperties(actual,expected);
        end

        function identicalRetentionArgumentsChooseNonfrozenBranch(testCase)
            p=testCase.fixture(); T=[-20;-2;-.01];
            h=p.L/(p.g*p.Tm)*T; [v,schema]=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyEqual(actual(:,2),zeros(3,1));
            testCase.verifyProperties(actual,expected);
        end

        function positiveHeadAddsTheExistingElasticInventory(testCase)
            p=testCase.fixture(); T=[0;4;110]; h=7;
            [v,schema]=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyEqual(actual(:,8),repmat(7*p.Ss,3,1));
            testCase.verifyProperties(actual,expected);
        end

        function scalarTemperatureExpansionMatchesPublicLaw(testCase)
            p=testCase.fixture(); T=-1; h=[-100;-1;0;5];
            [v,schema]=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifySize(actual,[4,12]);
            testCase.verifyProperties(actual,expected);
        end

        function missingOptionalSettingsRetainOriginalDefaults(testCase)
            p=rmfield(testCase.fixture(),'soilViscosity');
            [v,schema]=frozen_soil_jit_parameters(p); [T,h]=testCase.grid();
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyEqual(v(20),0); testCase.verifyEqual(v(22),0);
            testCase.verifyProperties(actual,expected);
        end

        function viscosityDisabledDoesNotRestrictTemperature(testCase)
            p=testCase.fixture(); p.soilViscosity.enabled=false;
            T=[-150;150]; h=[-2;2]; [v,schema]=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyProperties(actual,expected);
        end

        function zeroHydraulicConductivityStaysExactlyZero(testCase)
            p=testCase.fixture(); p.Ksat=0;
            [v,~]=frozen_soil_jit_parameters(p); [T,h]=testCase.grid();
            actual=frozen_soil_properties_jit(T,h,v);
            testCase.verifyEqual(actual(:,12),zeros(numel(T),1));
        end

        function iceImpedanceIsNotFlooredOrRemoved(testCase)
            p=testCase.fixture(); p.iceImpedance=200;
            [v,schema]=frozen_soil_jit_parameters(p); T=[-20;-2;4]; h=zeros(3,1);
            actual=frozen_soil_properties_jit(T,h,v);
            expected=testCase.reference(T,h,p,schema);
            testCase.verifyProperties(actual,expected);
            testCase.verifyLessThan(actual(1,12),1e-100);
        end

        function emptyPointArraysProduceAnEmptyPropertyMatrix(testCase)
            p=testCase.fixture(); v=frozen_soil_jit_parameters(p);
            actual=frozen_soil_properties_jit(zeros(0,1),zeros(0,1),v);
            testCase.verifySize(actual,[0,12]);
        end

        function schemaDescribesTheFixedScientificInterface(testCase)
            p=testCase.fixture(); [v,schema]=frozen_soil_jit_parameters(p);
            testCase.verifySize(v,[22,1]);
            testCase.verifyEqual(schema.propertyNames, ...
                {'sl','si','sg','thetaW','thetaI','thetaA','total', ...
                'elasticWater','mass','energy','k','K'});
            testCase.verifyEqual(v(20:22),[1;20;0]);
        end
    end

    methods (Access=private)
        function verifyProperties(testCase,actual,expected)
            testCase.verifySize(actual,size(expected));
            testCase.verifyEqual(actual(:,1:11),expected(:,1:11), ...
                'RelTol',4e-13,'AbsTol',2e-14);
            % A small absolute tolerance on K would conceal errors in frozen K.
            testCase.verifyEqual(actual(:,12),expected(:,12), ...
                'RelTol',4e-13,'AbsTol',1e-300);
        end
    end

    methods (Static,Access=private)
        function p=fixture()
            p=struct('rho',1000,'g',9.81,'L',3.34e5,'Tm',273.15, ...
                'porosity',.35,'thetaR',.05,'vgAlpha',.65,'vgN',1.67, ...
                'Ksat',9.81e-8,'Ss',1e-5,'iceImpedance',3.2, ...
                'cw',4180,'ci',2100,'rhoS',2650,'cs',800, ...
                'kw',.57,'ki',2.2,'ks',2.3,'ka',.025, ...
                'soilViscosity',struct('enabled',true,'referenceTemperature_C',20));
        end

        function [T,h]=grid()
            [T,h]=ndgrid([-20;-10;-1;-.01;-eps;0;eps;4;110], ...
                [-1e200;-1e8;-100;-1;-.001;-eps;0;eps;7]);
            T=T(:); h=h(:);
        end

        function Y=reference(T,h,p,schema)
            a=frozen_soil_constitutive(T,h,p,false);
            Y=zeros(numel(a.mass),12);
            for j=1:12, Y(:,j)=a.(schema.propertyNames{j})(:); end
        end
    end
end
