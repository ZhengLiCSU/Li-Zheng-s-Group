classdef TestJointGeometryParameterization < matlab.unittest.TestCase
    methods (Test)
        function inventoryInvertsGeneratedConductance(testCase)
            geometry = TestJointGeometryParameterization.fixture();
            freeze = joint_freeze_defaults();
            [conductance,reference] = joint_conductance_from_geometry(geometry,freeze);
            [~,~,inventory] = joint_channel_inventory(5,freeze, ...
                struct('enabled',true,'tortuosity',geometry.tortuosity),conductance, ...
                geometry.liningThickness_m,1000,3.35e5,4180,2100,2400*900,9.81);
            testCase.verifyEqual(inventory.waterVolume_m2,reference.waterVolume_m2,'RelTol',1e-13);
        end
        function radiusGrowthAtFixedVolumeChangesConductanceQuadratically(testCase)
            geometry = TestJointGeometryParameterization.fixture();
            freeze = joint_freeze_defaults();
            [first,a] = joint_conductance_from_geometry(geometry,freeze);
            freeze.flowWeightedRadius_m = 2*freeze.flowWeightedRadius_m;
            [second,b] = joint_conductance_from_geometry(geometry,freeze);
            testCase.verifyEqual(second,4*first,'RelTol',1e-14);
            testCase.verifyEqual(a.waterVolume_m2,b.waterVolume_m2,'AbsTol',0);
        end
        function longerPathsAtFixedVolumeDecreaseConductance(testCase)
            geometry = TestJointGeometryParameterization.fixture();
            first = joint_conductance_from_geometry(geometry,joint_freeze_defaults());
            geometry.tortuosity = 2*geometry.tortuosity;
            second = joint_conductance_from_geometry(geometry,joint_freeze_defaults());
            testCase.verifyEqual(second,first/4,'RelTol',1e-14);
        end
        function sealedFractionProducesZeroConductance(testCase)
            geometry = TestJointGeometryParameterization.fixture();
            geometry.equivalentVoidFraction = [0;0.01];
            [conductance,result] = joint_conductance_from_geometry(geometry,joint_freeze_defaults());
            testCase.verifyEqual(conductance(1),0,'AbsTol',0);
            testCase.verifyEqual(result.waterVolume_m2(1),0,'AbsTol',0);
        end
        function excessiveFractionIsRejected(testCase)
            geometry = TestJointGeometryParameterization.fixture();
            geometry.equivalentVoidFraction = 1.01;
            testCase.verifyError(@()joint_conductance_from_geometry(geometry,joint_freeze_defaults()), ...
                'JointGeometry:VoidFraction');
        end
        function negativeFractionIsRejected(testCase)
            geometry = TestJointGeometryParameterization.fixture();
            geometry.equivalentVoidFraction = -0.01;
            testCase.verifyError(@()joint_conductance_from_geometry(geometry,joint_freeze_defaults()), ...
                'JointGeometry:VoidFraction');
        end
    end
    methods (Static, Access=private)
        function geometry = fixture()
            geometry = struct('bandWidth_m',0.02,'liningThickness_m',0.45, ...
                'equivalentVoidFraction',0.01,'tortuosity',1.3, ...
                'density_kg_m3',1000,'gravity_m_s2',9.81);
        end
    end
end
