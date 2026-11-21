classdef TestSegmentContactFEM < matlab.unittest.TestCase
    %TESTSEGMENTCONTACTFEM Ring contact/bolt and interface algebra tests.
    % Run results=runtests('TestSegmentContactFEM.m'); assertSuccess(results).
    % These checks establish numerical mechanics, not gasket calibration.
    properties (Access=private)
        Source
        Mesh
        Cache
        ReferenceTemperature
    end
    methods (TestClassSetup)
        function buildGeometryAndCache(testCase)
            root=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.Source=load(fullfile(root,'example_mesh_repaired_interface.mat'));
            layout=hsr_eight_segment_layout();
            testCase.Mesh=build_segment_mechanical_mesh(testCase.Source,layout.jointAnglesRad);
            testCase.Cache=solve_segment_contact_fem('assemble',testCase.Mesh,struct());
            testCase.ReferenceTemperature=20*ones(size(testCase.Source.pAll,1),1);
        end
    end
    methods (Test)
        function zeroReferenceLoadHasZeroEquilibrium(testCase)
            out=solve_segment_contact_fem(testCase.Mesh,testCase.ReferenceTemperature, ...
                struct('boltPreload_N_per_m',0));
            testCase.verifyEqual(out.U_m,zeros(size(out.U_m)),'AbsTol',0);
            testCase.verifyEqual(out.jointGap_m,zeros(size(out.jointGap_m)),'AbsTol',0);
            testCase.verifyEqual(out.jointContactPressure_Pa,zeros(size(out.jointGap_m)),'AbsTol',0);
            testCase.verifyEqual(out.maximumFreeResidual_N_per_m,0,'AbsTol',0);
        end

        function eachElasticBlockHasThreeFreeRigidMotions(testCase)
            basis=TestSegmentContactFEM.rigidBasis(testCase.Mesh);
            K=testCase.Cache.Kelastic;
            error=norm(K*basis,'fro')/max(1,norm(K,'fro')*norm(basis,'fro'));
            testCase.verifyLessThan(error,2e-12);
            testCase.verifyEqual(rank(full(basis)),24);
        end

        function uniformPrescribedCompressionActivatesContact(testCase)
            M=testCase.Mesh;
            support=-1e-4*(M.p-M.center)/max(M.sourceRingRadii_m);
            out=solve_segment_contact_fem(M,testCase.ReferenceTemperature, ...
                struct('boltPreload_N_per_m',0,'foundationDisplacement_m',support));
            testCase.verifyGreaterThan(max(out.jointContactPressureGauss_Pa),0);
            testCase.verifyLessThan(min(out.jointGapGauss_m),0);
            testCase.verifyGreaterThanOrEqual(out.jointContactPressure_Pa,0);
            testCase.verifyLessThan(out.normalizedFreeResidual,2e-9);
        end

        function openedJointNeverCarriesTensileContact(testCase)
            M=testCase.Mesh;
            opened=TestSegmentContactFEM.radialTranslations(M,1e-3);
            out=solve_segment_contact_fem(M,testCase.ReferenceTemperature, ...
                struct('solveEquilibrium',false,'initialDisplacement_m',opened, ...
                'boltPreload_N_per_m',0,'boltNormalStiffness_N_per_m2',0));
            testCase.verifyGreaterThan(out.jointGapGauss_m,0);
            testCase.verifyEqual(out.jointContactPressureGauss_Pa, ...
                zeros(size(out.jointContactPressureGauss_Pa)),'AbsTol',0);
        end

        function boltPreloadBalancesContactAndSupport(testCase)
            out=solve_segment_contact_fem(testCase.Mesh,testCase.ReferenceTemperature,struct());
            testCase.verifyGreaterThan(min(out.boltForce_N_per_m),0);
            testCase.verifyGreaterThan(max(out.jointContactPressureGauss_Pa),0);
            testCase.verifyLessThan(out.normalizedFreeResidual,2e-9);
            ringInterfaceForce=-testCase.Cache.KringRing*out.uVector_m;
            soilForce=out.soilReactionForce_global_N_per_m;
            resultant=sum(reshape(ringInterfaceForce,2,[])',1)+sum(soilForce,1);
            testCase.verifyEqual(resultant,[0 0],'AbsTol',2e-6);
        end

        function thermalLoadProducesActualGapAndFaceRotation(testCase)
            temperature=testCase.ReferenceTemperature+5+.4*testCase.Source.pAll(:,2);
            out=solve_segment_contact_fem(testCase.Mesh,temperature, ...
                struct('boltPreload_N_per_m',0));
            expected=TestSegmentContactFEM.gapsFromDisplacement(testCase.Mesh,out.U_m);
            testCase.verifyGreaterThan(norm(out.U_m,'fro'),1e-8);
            testCase.verifyGreaterThan(max(abs(out.jointGap_m),[],'all'),1e-12);
            testCase.verifyGreaterThan(max(abs(out.relativeFaceRotation_rad)),1e-12);
            testCase.verifyEqual(out.jointGap_m,expected,'AbsTol',1e-13);
        end

        function fullyConstrainedThermalStressUsesThreeBulkModulus(testCase)
            M=testCase.Mesh;
            out=solve_segment_contact_fem(M,testCase.ReferenceTemperature+1, ...
                struct('boltPreload_N_per_m',0,'fixedDofs',(1:2*size(M.p,1))'));
            coefficient=out.parameters.E_Pa/(1-2*out.parameters.poissonRatio)* ...
                out.parameters.thermalExpansion_per_C;
            expected=-coefficient*ones(size(M.q,1),4);
            testCase.verifyEqual(out.stressGauss_Pa(:,:,1),expected,'AbsTol',2e-5);
            testCase.verifyEqual(out.stressGauss_Pa(:,:,2),expected,'AbsTol',2e-5);
            testCase.verifyEqual(out.stressGauss_Pa(:,:,4),expected,'AbsTol',2e-5);
            testCase.verifyEqual(out.stressGauss_Pa(:,:,3),zeros(size(expected)),'AbsTol',1e-10);
        end

        function cachedContactAndSoilTangentMatchesIndependentDifference(testCase)
            [state,direction]=TestSegmentContactFEM.mixedState(testCase.Mesh,testCase.Source);
            out=solve_segment_contact_fem('evaluate',testCase.Cache,state);
            exact=[out.ringTangent_N_per_m2,out.ringSoilTangent_N_per_m2; ...
                out.soilRingTangent_N_per_m2,out.soilSoilTangent_N_per_m2]*direction;
            observed=TestSegmentContactFEM.directionalDifference(testCase.Cache,state,direction,1e-8);
            testCase.verifyTrue(any(out.jointGapGauss_m<0)&&any(out.jointGapGauss_m>0));
            testCase.verifyLessThan(norm(observed-exact)/max(1,norm(exact)),2e-7);
        end

        function cachedThermalLoadMapMatchesIndependentDifference(testCase)
            [state,~]=TestSegmentContactFEM.mixedState(testCase.Mesh,testCase.Source);
            stream=RandStream('mt19937ar','Seed',771);
            direction=randn(stream,numel(state.temperatureGlobal_C),1);
            direction=direction/norm(direction);
            plus=state; minus=state;
            plus.temperatureGlobal_C=state.temperatureGlobal_C+1e-3*direction;
            minus.temperatureGlobal_C=state.temperatureGlobal_C-1e-3*direction;
            a=solve_segment_contact_fem('evaluate',testCase.Cache,plus);
            b=solve_segment_contact_fem('evaluate',testCase.Cache,minus);
            finite=(a.ringResidual_N_per_m-b.ringResidual_N_per_m)/(2e-3);
            exact=-testCase.Cache.thermalLoadMap*direction;
            testCase.verifyLessThan(norm(finite-exact)/max(1,norm(exact)),2e-7);
        end

        function cachedInterfaceHasSymmetricActionReactionBlocks(testCase)
            C=testCase.Cache;
            testCase.verifyEqual(C.KringSoil,-C.Wring'*C.Dinterface*C.Wsoil,'AbsTol',1e-8);
            testCase.verifyEqual(C.KsoilSoil,C.Wsoil'*C.Dinterface*C.Wsoil,'AbsTol',1e-8);
            uniformRing=repmat([.001 -.002],size(testCase.Mesh.p,1),1);
            uniformSoil=repmat([.001 -.002],size(testCase.Source.pAll,1),1);
            state=struct('temperatureGlobal_C',testCase.ReferenceTemperature, ...
                'ringDisplacement_m',uniformRing,'soilDisplacementGlobal_m',uniformSoil);
            out=solve_segment_contact_fem('evaluate',C,state);
            testCase.verifyEqual(out.soilResidual_N_per_m,zeros(size(out.soilResidual_N_per_m)), ...
                'AbsTol',2e-7);
        end

        function cachedAndFreshEvaluationAgree(testCase)
            [state,~]=TestSegmentContactFEM.mixedState(testCase.Mesh,testCase.Source);
            state.computeStress=true;
            cached=solve_segment_contact_fem('evaluate',testCase.Cache,state);
            fresh=solve_segment_contact_fem(testCase.Mesh,state.temperatureGlobal_C, ...
                struct('solveEquilibrium',false,'initialDisplacement_m',state.ringDisplacement_m, ...
                'foundationDisplacementGlobal_m',state.soilDisplacementGlobal_m));
            testCase.verifyEqual(cached.ringResidual_N_per_m,fresh.residual_N_per_m,'AbsTol',2e-6);
            testCase.verifyEqual(cached.soilResidual_N_per_m,fresh.soilInterfaceResidual_N_per_m, ...
                'AbsTol',2e-6);
            testCase.verifyEqual(cached.stressGauss_Pa,fresh.stressGauss_Pa,'AbsTol',1e-8);
        end

        function unsupportedRigidModesAreRejected(testCase)
            settings=struct('foundationNormalStiffness_Pa_per_m',0, ...
                'foundationTangentialStiffness_Pa_per_m',0,'boltPreload_N_per_m',0);
            testCase.verifyError(@() solve_segment_contact_fem(testCase.Mesh, ...
                testCase.ReferenceTemperature,settings),'SegmentContactFEM:RigidModes');
        end

        function cachedEvaluationPreservesMechanicalSupportTarget(testCase)
            [state,~]=TestSegmentContactFEM.mixedState(testCase.Mesh,testCase.Source);
            state=rmfield(state,'soilDisplacementGlobal_m');
            support=-1e-4*(testCase.Mesh.p-testCase.Mesh.center)/max(testCase.Mesh.sourceRingRadii_m);
            settings=struct('foundationDisplacement_m',support);
            cache=solve_segment_contact_fem('assemble',testCase.Mesh,settings);
            cached=solve_segment_contact_fem('evaluate',cache,state);
            settings.solveEquilibrium=false;
            settings.initialDisplacement_m=state.ringDisplacement_m;
            fresh=solve_segment_contact_fem(testCase.Mesh,state.temperatureGlobal_C,settings);
            testCase.verifyEqual(cached.ringResidual_N_per_m,fresh.residual_N_per_m,'AbsTol',2e-6);
            testCase.verifyEqual(cached.soilResidual_N_per_m,fresh.soilInterfaceResidual_N_per_m, ...
                'AbsTol',2e-6);
        end

        function cachedEvaluationPreservesGlobalSoilSupportTarget(testCase)
            [state,~]=TestSegmentContactFEM.mixedState(testCase.Mesh,testCase.Source);
            support=state.soilDisplacementGlobal_m;
            state=rmfield(state,'soilDisplacementGlobal_m');
            settings=struct('foundationDisplacementGlobal_m',support);
            cache=solve_segment_contact_fem('assemble',testCase.Mesh,settings);
            cached=solve_segment_contact_fem('evaluate',cache,state);
            settings.solveEquilibrium=false;
            settings.initialDisplacement_m=state.ringDisplacement_m;
            fresh=solve_segment_contact_fem(testCase.Mesh,state.temperatureGlobal_C,settings);
            testCase.verifyEqual(cached.ringResidual_N_per_m,fresh.residual_N_per_m,'AbsTol',2e-6);
            testCase.verifyEqual(cached.soilResidual_N_per_m,fresh.soilInterfaceResidual_N_per_m, ...
                'AbsTol',2e-6);
        end
    end

    methods (Static,Access=private)
        function basis=rigidBasis(M)
            N=size(M.p,1); basis=sparse(2*N,24);
            for segment=1:8
                nodes=find(M.nodeSegmentID==segment); dofX=2*nodes-1; dofY=2*nodes;
                column=3*(segment-1); centre=mean(M.p(nodes,:),1);
                basis(dofX,column+1)=1; basis(dofY,column+2)=1;
                basis(dofX,column+3)=-(M.p(nodes,2)-centre(2));
                basis(dofY,column+3)=M.p(nodes,1)-centre(1);
            end
        end

        function u=radialTranslations(M,magnitude)
            u=zeros(size(M.p));
            for segment=1:8
                nodes=find(M.nodeSegmentID==segment);
                direction=mean(M.p(nodes,:),1)-M.center; direction=direction/norm(direction);
                u(nodes,:)=repmat(magnitude*direction,numel(nodes),1);
            end
        end

        function gap=gapsFromDisplacement(M,U)
            gap=zeros(size(M.jointNodePairs,1),8);
            for j=1:8
                jump=U(M.jointNodePairs(:,2,j),:)-U(M.jointNodePairs(:,1,j),:);
                gap(:,j)=jump*[-sin(M.jointAnglesRad(j));cos(M.jointAnglesRad(j))];
            end
        end

        function [state,direction]=mixedState(M,S)
            stream=RandStream('mt19937ar','Seed',121);
            ring=1e-4*randn(stream,size(M.p)); soil=1e-5*randn(stream,size(S.pAll));
            state=struct('temperatureGlobal_C',20+.1*S.pAll(:,2), ...
                'ringDisplacement_m',ring,'soilDisplacementGlobal_m',soil);
            direction=randn(stream,2*size(M.p,1)+2*size(S.pAll,1),1);
            direction=direction/norm(direction);
        end

        function difference=directionalDifference(cache,state,direction,step)
            nr=cache.mechanicalRingNodeCount;
            ring=reshape(state.ringDisplacement_m',[],1);
            soil=reshape(state.soilDisplacementGlobal_m',[],1);
            plus=state; minus=state;
            plus.ringDisplacement_m=ring+step*direction(1:2*nr);
            minus.ringDisplacement_m=ring-step*direction(1:2*nr);
            plus.soilDisplacementGlobal_m=soil+step*direction(2*nr+1:end);
            minus.soilDisplacementGlobal_m=soil-step*direction(2*nr+1:end);
            a=solve_segment_contact_fem('evaluate',cache,plus);
            b=solve_segment_contact_fem('evaluate',cache,minus);
            difference=([a.ringResidual_N_per_m;a.soilResidual_N_per_m]- ...
                [b.ringResidual_N_per_m;b.soilResidual_N_per_m])/(2*step);
        end
    end
end
