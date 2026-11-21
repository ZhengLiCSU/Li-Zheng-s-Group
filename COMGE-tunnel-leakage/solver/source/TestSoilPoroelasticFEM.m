classdef TestSoilPoroelasticFEM < matlab.unittest.TestCase
    % Mechanics patches independent of tunnel contact and displacement BCs.
    methods (Test)
        function rigidTranslationsAndRotationHaveNoElasticForce(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            a=soil_poroelastic_fem(xy,tri,5,3,5,s); n=size(xy,1);
            tx=reshape([ones(1,n);zeros(1,n)],[],1);
            tz=reshape([zeros(1,n);ones(1,n)],[],1);
            rotation=reshape([-xy(:,2)';xy(:,1)'],[],1);
            testCase.verifyEqual(a.K,a.K','AbsTol',2e-8);
            scale=max(1,norm(a.K,inf));
            testCase.verifyLessThan(norm(a.K*[tx tz rotation],inf)/scale,2e-14);
            testCase.verifyEqual(a.V*[tx tz rotation],zeros(size(tri,1),3),'AbsTol',1e-14);
        end

        function affineDilationAndSCVProjectionAreExact(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            a=soil_poroelastic_fem(xy,tri,5,3,5,s);
            u=reshape([.013*xy(:,1)';-.004*xy(:,2)'],[],1);
            testCase.verifyEqual(a.V*u,.009*ones(size(tri,1),1),'AbsTol',1e-15);
            testCase.verifyEqual(a.nodalVolumetricStrainMatrix*u, ...
                .009*ones(size(xy,1),1),'AbsTol',2e-15);
            testCase.verifyEqual(a.nodeProjection*ones(size(tri,1),1), ...
                ones(size(xy,1),1),'AbsTol',2e-15);
            testCase.verifyEqual(sum(a.nodeVolumes),1,'AbsTol',2e-15);
        end

        function affineStressMatchesThreeDimensionalPlaneStrainLaw(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            a=soil_poroelastic_fem(xy,tri,5,3,5,s);
            ux=.012*xy(:,1)-.003*xy(:,2)+.8;
            uz=.004*xy(:,1)-.006*xy(:,2)-.3;
            u=reshape([ux';uz'],[],1); strain=[.012;-.006;.001];
            G=s.E/(2*(1+s.nu)); lambda=s.E*s.nu/((1+s.nu)*(1-2*s.nu));
            expected=[2*G*strain(1)+lambda*(strain(1)+strain(2)); ...
                2*G*strain(2)+lambda*(strain(1)+strain(2));G*strain(3)];
            stored=0;
            for e=1:size(tri,1)
                ids=tri(e,:); dof=reshape([2*ids-1;2*ids],[],1);
                actualStrain=a.B(:,:,e)*u(dof);
                testCase.verifyEqual(actualStrain,strain,'AbsTol',2e-15);
                testCase.verifyEqual(a.D(:,:,e)*actualStrain,expected,'AbsTol',1e-7);
                stored=stored+.5*a.elementAreas(e)*strain'*expected;
            end
            testCase.verifyEqual(.5*u'*a.K*u,stored,'AbsTol',1e-7,'RelTol',1e-11);
        end

        function thermalWorkIncludesConstrainedOutOfPlaneEigenstress(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            dT=7; a=soil_poroelastic_fem(xy,tri,5+dT,3,5,s);
            u=reshape([.01*xy(:,1)';.02*xy(:,2)'],[],1);
            beta=s.E*s.alphaThermal/(1-2*s.nu);
            expectedWork=beta*dT*.03; % unit area, constant div(u)
            testCase.verifyEqual(a.thermalEigenstressCoefficient_Pa_K, ...
                beta*ones(size(tri,1),1),'AbsTol',1e-10);
            testCase.verifyEqual(u'*a.fthermal,expectedWork,'AbsTol',1e-9,'RelTol',1e-13);
            testCase.verifyGreaterThan(u'*a.fthermal,0);
        end

        function poreLoadHasCorrectBoundaryWorkAndSign(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture(); s.phasePressure='liquid';
            h=2+xy(:,2); a=soil_poroelastic_fem(xy,tri,5,h,5,s);
            u=reshape([.01*xy(:,1)';.02*xy(:,2)'],[],1);
            expected=s.biot*s.phaseParameters.rho*s.phaseParameters.g*2.5*.03;
            testCase.verifyEqual(u'*a.fpore,expected,'AbsTol',1e-10,'RelTol',1e-13);
            testCase.verifyGreaterThan(u'*a.fpore,0);
            net=[sum(a.fpore(1:2:end));sum(a.fpore(2:2:end))];
            testCase.verifyEqual(net,[0;0],'AbsTol',1e-10);
        end

        function hydrostaticGradientAndBodyForceBalanceWithTractions(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture(); s.phasePressure='liquid';
            s.rhoBulk=s.biot*s.phaseParameters.rho;
            h=1+xy(:,2); a=soil_poroelastic_fem(xy,tri,5,h,5,s);
            force=zeros(2*size(xy,1),1); edges=[1 2;2 4;4 3;3 1];
            for e=1:size(edges,1)
                ids=edges(e,:); tangent=xy(ids(2),:)-xy(ids(1),:); L=norm(tangent);
                normal=[tangent(2);-tangent(1)]/L;
                pressure=s.phaseParameters.rho*s.phaseParameters.g*h(ids);
                traction=-s.biot*normal*pressure';
                load=L/6*traction*[2 1;1 2];
                dof=reshape([2*ids-1;2*ids],[],1);
                force(dof)=force(dof)+load(:);
            end
            testCase.verifyEqual(a.fbody+a.fpore+force,zeros(size(force)), ...
                'AbsTol',2e-10);
            testCase.verifyEqual(sum(a.fbody(1:2:end)),0,'AbsTol',0);
            testCase.verifyEqual(sum(a.fbody(2:2:end)),s.rhoBulk*s.gravity(2), ...
                'AbsTol',1e-10);
        end

        function saturatedClapeyronPressureAndHeadTangentAreExact(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture(); T=-.4; h=3;
            a=soil_poroelastic_fem(xy,tri,T,h,5,s); q=a.quadrature;
            p=s.phaseParameters; jump=-p.rho*p.L/p.Tm*T;
            testCase.verifyEqual(q.piCandidate_Pa-q.pw_Pa, ...
                jump*ones(size(q.T)),'AbsTol',1e-9);
            testCase.verifyEqual(q.pbar_Pa,q.pw_Pa+q.si*jump, ...
                'AbsTol',1e-9,'RelTol',1e-14);
            testCase.verifyEqual(q.dpbarH_Pa_m,p.rho*p.g*ones(size(q.T)), ...
                'AbsTol',1e-9);
            testCase.verifyLessThan(q.dpbarT_Pa_K,zeros(size(q.T)));
            epsT=2e-7;
            plus=soil_poroelastic_fem(xy,tri,T+epsT,h,5,s);
            minus=soil_poroelastic_fem(xy,tri,T-epsT,h,5,s);
            fd=(plus.quadrature.pbar_Pa-minus.quadrature.pbar_Pa)/(2*epsT);
            testCase.verifyEqual(q.dpbarT_Pa_K,fd,'AbsTol',.02,'RelTol',2e-6);
        end

        function completeLoadTangentsMatchDirectionalDifferences(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            T=[-1.2;-.8;-1;-.9;-1.1]; h=[-2;-1;3;-4;.2];
            directionT=[-.2;.3;.4;-.1;.1]; directionH=[.4;-.2;.1;.3;-.1]; epsilon=2e-7;
            a=soil_poroelastic_fem(xy,tri,T,h,5,s);
            tp=soil_poroelastic_fem(xy,tri,T+epsilon*directionT,h,5,s);
            tm=soil_poroelastic_fem(xy,tri,T-epsilon*directionT,h,5,s);
            hp=soil_poroelastic_fem(xy,tri,T,h+epsilon*directionH,5,s);
            hm=soil_poroelastic_fem(xy,tri,T,h-epsilon*directionH,5,s);
            fdt=(tp.fthermal+tp.fpore-tm.fthermal-tm.fpore)/(2*epsilon);
            fdh=(hp.fthermal+hp.fpore-hm.fthermal-hm.fpore)/(2*epsilon);
            testCase.verifyEqual(a.dloadT*directionT,fdt,'AbsTol',.03,'RelTol',3e-6);
            testCase.verifyEqual(a.dloadH*directionH,fdh,'AbsTol',.03,'RelTol',3e-6);
        end

        function cachedLoadMapsAndQuadratureReproduceLoads(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            T=[-1;-.5;-.9;-.7;-.8]; h=[-2;3;-1;4;.1];
            a=soil_poroelastic_fem(xy,tri,T,h,5,s);
            testCase.verifyEqual(a.quadrature.P*T,a.quadrature.T(:),'AbsTol',2e-16);
            testCase.verifyEqual(a.quadrature.P*h,a.quadrature.h(:),'AbsTol',1e-15);
            testCase.verifyEqual(a.A_temperature*(a.elementMeanNodal*(T-5)), ...
                a.fthermal,'AbsTol',1e-9,'RelTol',2e-14);
            testCase.verifyEqual(a.A_pressure*(a.quadrature.elementAverage* ...
                a.quadrature.pbar_Pa(:)),a.fpore,'AbsTol',1e-9,'RelTol',2e-14);
        end

        function triangleOrientationDoesNotChangePhysicalOperators(testCase)
            [xy,tri,s]=TestSoilPoroelasticFEM.fixture();
            a=soil_poroelastic_fem(xy,tri,-.4,3,5,s);
            b=soil_poroelastic_fem(xy,tri(:,[1 3 2]),-.4,3,5,s);
            testCase.verifyEqual(a.K,b.K,'AbsTol',2e-8,'RelTol',2e-14);
            testCase.verifyEqual(a.V,b.V,'AbsTol',1e-15);
            testCase.verifyEqual(a.fpore,b.fpore,'AbsTol',1e-9);
            testCase.verifyEqual(a.fthermal,b.fthermal,'AbsTol',1e-9);
        end
    end
    methods (Static,Access=private)
        function [xy,tri,s]=fixture()
            xy=[0 0;1 0;0 1;1 1;.5 .5];
            tri=[1 2 5;2 4 5;4 3 5;3 1 5];
            p=struct('rho',1000,'g',9.81,'L',3.34e5,'Tm',273.15, ...
                'porosity',.35,'thetaR',.05,'vgAlpha',.65,'vgN',1.67, ...
                'Ksat',9.81e-8,'Ss',1e-5,'iceImpedance',0, ...
                'cw',4180,'ci',2100,'rhoS',2650,'cs',800, ...
                'kw',.57,'ki',2.2,'ks',2.3,'ka',.025, ...
                'soilViscosity',struct('enabled',false,'referenceTemperature_C',20));
            s=struct('E',80e6,'nu',.35,'alphaThermal',1.2e-5,'biot',.65, ...
                'rhoBulk',2000,'gravity',[0;9.81],'phasePressure','bishopIce', ...
                'phaseParameters',p);
        end
    end
end
