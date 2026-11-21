classdef TestJointDrainQuadrature < matlab.unittest.TestCase
    %TESTJOINTDRAINQUADRATURE Independent public-interface Robin checks.
    % The diamond interface gives exact P1 edge integrals. It does not
    % require a tunnel solver, existing result files or a MATLAB desktop.
    properties (TestParameter)
        supportWidth=struct('narrow',.02,'medium',.4,'wholeEdge',pi/2);
    end

    methods (Test)
        function uniformHeadPreservesIntegratedConductance(testCase,supportWidth)
            [S,L,B,p]=TestJointDrainQuadrature.fixture(pi/4,supportWidth,false);
            Q=build_joint_drain_quadrature(S,(1:4)',L,[0 pi/2 pi 3*pi/2]);
            B.jointQuadrature=Q;
            [q,ledger]=evaluate_joint_drain(4*ones(4,1),5*ones(4,1),B,p,1);
            testCase.verifyEqual(sum(q),8,'AbsTol',1e-12);
            testCase.verifyEqual(ledger.byJoint,8,'AbsTol',1e-12);
            testCase.verifyEqual(sum(ledger.bySegment),8,'AbsTol',1e-12);
            testCase.verifyEqual(full(sum(Q.phi,2)),ones(size(Q.phi,1),1),'AbsTol',1e-13);
            testCase.verifyEqual(sum(Q.weight),1,'AbsTol',1e-13);
            testCase.verifyEqual(Q.phi*S.pAll,Q.sampleXY,'AbsTol',1e-13);
            testCase.verifyGreaterThan(Q.supportLength,0);
        end

        function linearHeadHasConsistentEdgeIntegral(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,pi/2,false);
            [q,ledger,DH]=evaluate_joint_drain([2;6;2;2],5*ones(4,1),B,p,1);
            expected=zeros(4); expected(1:2,1:2)=[2 1;1 2]/3;
            testCase.verifyEqual(q,[10/3;14/3;0;0],'AbsTol',2e-12);
            testCase.verifyEqual(full(DH),expected,'AbsTol',1e-12);
            testCase.verifyEqual(full(ledger.headSlope),[1 1 0 0],'AbsTol',1e-12);
            testCase.verifyEqual(ledger.drivingHead,4,'AbsTol',1e-12);
            testCase.verifyEqual(ledger.bySegment,[8;0;0;0],'AbsTol',1e-12);
            testCase.verifyEqual(B.jointQuadrature.supportLength,sqrt(2),'AbsTol',1e-12);
        end

        function headDerivativeMatchesFiniteDifference(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,pi/2,false);
            h=[2;6;2;2]; T=5*ones(4,1); v=[.2;-.3;.1;-.1]; epsilon=1e-6;
            [~,~,DH]=evaluate_joint_drain(h,T,B,p,1);
            qp=evaluate_joint_drain(h+epsilon*v,T,B,p,1);
            qm=evaluate_joint_drain(h-epsilon*v,T,B,p,1);
            testCase.verifyEqual(DH*v,(qp-qm)/(2*epsilon),'AbsTol',1e-9);
        end

        function liningTemperatureDerivativeIsNonlocal(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,pi/2,true);
            h=[2;6;2;2]; T=[2;4;1;3;-.06;-.04;4;5];
            v=[.15;-.1;.05;.1;-.003;.004;.02;-.03]; epsilon=1e-5;
            [~,ledger,~,DT]=evaluate_joint_drain(h,T,B,p,1);
            qp=evaluate_joint_drain(h,T+epsilon*v,B,p,1);
            qm=evaluate_joint_drain(h,T-epsilon*v,B,p,1);
            testCase.verifyEqual(DT*v,(qp-qm)/(2*epsilon),'AbsTol',1e-8,'RelTol',1e-6);
            testCase.verifyGreaterThan(norm(DT(:,5:6),1),0);
            testCase.verifySize(ledger.dFactorT,[1 8]);
            testCase.verifyEqual(qp(3:4),zeros(2,1),'AbsTol',0);
        end

        function inletHeatHeadDerivativeMatchesFiniteDifference(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,pi/2,false);
            h=[2;6;2;2]; T=[2;5;7;3]; v=[.2;-.3;.1;-.1]; epsilon=1e-6;
            [~,ledger]=evaluate_joint_drain(h,T,B,p,1);
            [~,lp]=evaluate_joint_drain(h+epsilon*v,T,B,p,1);
            [~,lm]=evaluate_joint_drain(h-epsilon*v,T,B,p,1);
            fd=(lp.nodeFlow*T-lm.nodeFlow*T)/(2*epsilon);
            testCase.verifyEqual(ledger.entryHeatHeadDerivative*v,fd,'AbsTol',2e-9);
            testCase.verifyEqual(full(ledger.entryHeatHeadDerivative),[3 4 0 0],'AbsTol',1e-12);
        end

        function negativeDriveDisablesDrainAndDerivatives(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,pi/2,true);
            [q,ledger,DH,DT]=evaluate_joint_drain(-ones(4,1), ...
                [2;4;1;3;-.06;-.04;4;5],B,p,1);
            testCase.verifyEqual(q,zeros(4,1),'AbsTol',0);
            testCase.verifyEqual(ledger.byJoint,0,'AbsTol',0);
            testCase.verifyEqual(full(DH),zeros(4),'AbsTol',0);
            testCase.verifyEqual(full(DT),zeros(4,8),'AbsTol',0);
        end

        function bandAtSegmentBoundaryPartitionsWithoutLosingWater(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(0,pi/2,false);
            [q,ledger]=evaluate_joint_drain(4*ones(4,1),5*ones(4,1),B,p,1);
            testCase.verifyEqual(sum(q),8,'AbsTol',1e-12);
            testCase.verifyEqual(ledger.bySegment,[4;0;0;4],'AbsTol',1e-12);
            testCase.verifyEqual(sum(ledger.bySegment),ledger.byJoint,'AbsTol',1e-12);
        end

        function narrowBandApproachesSymmetricPositiveRankOne(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,1e-6,false);
            [~,~,DH]=evaluate_joint_drain(4*ones(4,1),5*ones(4,1),B,p,1);
            spectrum=eig(full(DH(1:2,1:2)));
            testCase.verifyEqual(full(DH),full(DH'),'AbsTol',1e-14);
            testCase.verifyGreaterThanOrEqual(min(eig(full(DH))),-1e-12);
            testCase.verifyEqual(full(DH(1:2,1:2)),.5*ones(2),'AbsTol',2e-10);
            testCase.verifyLessThan(spectrum(1)/spectrum(2),1e-10);
        end

        function staticHeadReferenceRetainsEntryAndAllowsStorageExcess(testCase)
            [~,~,B,p]=TestJointDrainQuadrature.fixture(pi/4,pi/2,false);
            p.leakage.jointEntryHead=1; B.soilDepth=[3;5;3;3];
            [q,ledger]=evaluate_joint_drain(4*ones(4,1),5*ones(4,1),B,p,1);
            testCase.verifyEqual(ledger.drivingHead,4,'AbsTol',1e-12);
            testCase.verifyEqual(ledger.referenceHead,3,'AbsTol',1e-12);
            testCase.verifyEqual(ledger.unrestrictedRate,4,'AbsTol',1e-12);
            testCase.verifyEqual(sum(q),6,'AbsTol',1e-12);
            testCase.verifyGreaterThan(sum(q),ledger.unrestrictedRate);
        end

        function circularSectorOrderingSurvivesRotation(testCase)
            [S,L,B,p]=TestJointDrainQuadrature.fixture(0,pi/2,false);
            B.jointQuadrature=build_joint_drain_quadrature(S,(1:4)',L, ...
                [3*pi/2 2*pi 5*pi/2 3*pi]);
            [~,ledger]=evaluate_joint_drain(4*ones(4,1),5*ones(4,1),B,p,1);
            testCase.verifyEqual(ledger.bySegment,[4;4;0;0],'AbsTol',1e-12);
        end
    end

    methods (Static,Access=private)
        function [S,L,B,p]=fixture(jointAngle,width,frozen)
            xy=[1 0;0 1;-1 0;0 -1]; edges=[1 2;2 3;3 4;4 1];
            table=zeros(4,7); table(:,2:3)=edges;
            S=struct('pAll',xy,'xc',0,'yc',0,'soilTunnelCoupleTable',table);
            L=struct('jointAnglesRad',jointAngle,'jointWidth_m',width, ...
                'jointConductance',2,'jointEntryHead',0,'interiorHead',0, ...
                'jointFreeze',joint_freeze_defaults(struct('enabled',frozen)));
            Q=build_joint_drain_quadrature(S,(1:4)',L,[0 pi/2 pi 3*pi/2]);
            B=struct('jointQuadrature',Q,'liningZDrop',zeros(4,1), ...
                'soilDepth',3*ones(4,1));
            if frozen
                B.soilSelect=sparse(1:4,1:4,1,4,8);
                B.jointThermalPaths=struct('nj',1,'nq',2, ...
                    'interpolation',sparse([1 2],[5 6],1,2,8));
            end
            p=struct('leakage',L,'rho',1000,'cw',4180);
        end
    end
end
