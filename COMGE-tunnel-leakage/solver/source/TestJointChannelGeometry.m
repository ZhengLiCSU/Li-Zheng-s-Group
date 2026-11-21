classdef TestJointChannelGeometry < matlab.unittest.TestCase
    % Reject nonphysical conductance/radius combinations before simulation.
    methods (Test)
        function validBandReportsEquivalentVoidFraction(testCase)
            result = audit_joint_channel_geometry(struct('waterVolume_m2',[0;0.001]),0.02,0.35);
            testCase.verifyEqual(result.equivalentVoidFraction,[0;1/7],'AbsTol',1e-15);
        end
        function excessInventoryIsRejected(testCase)
            testCase.verifyError(@()audit_joint_channel_geometry( ...
                struct('waterVolume_m2',0.008),0.02,0.35),'JointGeometry:ExcessVolume');
        end
        function dimensionsMustMatch(testCase)
            testCase.verifyError(@()audit_joint_channel_geometry( ...
                struct('waterVolume_m2',[0;0;0]),[0.01;0.02],0.35),'JointGeometry:Shape');
        end
        function closedJointUsesNoConnectedVolume(testCase)
            result = audit_joint_channel_geometry(struct('waterVolume_m2',zeros(8,1)),0.02,0.35);
            testCase.verifyEqual(result.equivalentVoidFraction,zeros(8,1));
        end
        function actualBaselineGeometryPasses(testCase)
            result = TestJointChannelGeometry.geometry(1e-6);
            testCase.verifyGreaterThan(result.equivalentVoidFraction(7),0);
            testCase.verifyLessThan(result.equivalentVoidFraction(7),0.01);
        end
        function excessivelySmallRadiusAtFixedConductanceFails(testCase)
            testCase.verifyError(@()TestJointChannelGeometry.geometry(0.02e-6), ...
                'JointGeometry:ExcessVolume');
        end
        function publicCoreChecksTheBandBeforeTheMeshVolume(testCase)
            model = joint_contact_thm_model(struct('mesh',make_thm_verification_mesh(), ...
                'durationDays',0.1,'saveOutput',false,'verbose',false,'plot',false, ...
                'leakage',struct('jointFreeze',struct('flowWeightedRadius_m',0.02e-6))));
            testCase.verifyError(@()TH_2D_segment_v16(model),'JointGeometry:ExcessVolume');
        end
    end
    methods (Static, Access=private)
        function result = geometry(radius)
            model = joint_supply_th_model(struct('mesh',make_thm_verification_mesh(), ...
                'leakage',struct('jointFreeze',struct('flowWeightedRadius_m',radius))));
            [~,~,inventory] = joint_channel_inventory(zeros(8,8),model.leakage.jointFreeze, ...
                model.leakage.jointStorage,model.leakage.jointConductance,0.35, ...
                1000,3.35e5,4180,2100,2400*900,9.81);
            result = audit_joint_channel_geometry(inventory,model.leakage.jointWidth_m,0.35);
        end
    end
end
