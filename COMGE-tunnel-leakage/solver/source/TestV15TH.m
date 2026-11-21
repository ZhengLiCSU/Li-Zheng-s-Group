classdef TestV15TH < matlab.unittest.TestCase
    methods (Test)
        function integratedCoreInvariants(testCase)
            report=TH_2D_segment_v15('selftest');
            testCase.verifyTrue(report.passed);
            testCase.verifyLessThan(report.enthalpyDatumPatchRelativeError,1e-10);
            testCase.verifyLessThan(report.realMeshJacobianRelativeError,1e-5);
            testCase.verifyLessThan(report.frozenJointJacobianRelativeError,1e-5);
            testCase.verifyLessThan(report.maxStepWaterDefect_m2,1e-7);
        end
        function checkpointPreservesStateAndRejectsOldVersion(testCase)
            report=TestV15TH.restartExercise();
            testCase.verifyEqual(report.timeDay,.3,'AbsTol',1e-13);
            testCase.verifyEqual(report.initialX,report.checkpointX,'AbsTol',0);
            testCase.verifyEqual(report.finalT,report.referenceT,'AbsTol',1e-8);
            testCase.verifyEqual(report.finalH,report.referenceH,'AbsTol',1e-8);
            testCase.verifyTrue(report.oldCheckpointRejected);
            testCase.verifyTrue(report.changedQuadratureRejected);
        end
    end
    methods (Static,Access=private)
        function report=restartExercise()
            folder=tempname(fileparts(mfilename('fullpath'))); mkdir(folder);
            mesh=struct('pAll',[0 0;1 0;0 1;1 1],'tSoil',[1 2 3;2 4 3]);
            o=struct('mesh',mesh,'durationDays',.3,'dtInitial',.05*86400, ...
                'dtMax',.05*86400,'outputEveryDays',.15,'topMin',1,'topMax',1, ...
                'hydroFar',false,'topSeepage',false,'tunnelDrain',false, ...
                'geothermalFlux',0,'advection',true,'saveOutput',false,'plot',false, ...
                'verbose',false,'mechanics',struct('enable',false), ...
                'initial',struct('surfaceT',2,'gradient',0), ...
                'solver',struct('predictor',false,'jacobianReuse',0));
            reference=TH_2D_segment_v15(o);
            o.durationDays=.15; o.saveOutput=true; o.outputPrefix=fullfile(folder,'partial');
            first=TH_2D_segment_v15(o); file=[o.outputPrefix '_checkpoint.mat'];
            d=load(file,'checkpoint');
            o.initial=struct('mode','checkpoint','checkpointFile',file);
            continuation=TH_2D_segment_v15(o);
            changed=o; changed.integration=struct('soilQuadratureOrder',2);
            changedRejected=false;
            try
                TH_2D_segment_v15(changed);
            catch errorInfo
                changedRejected=contains(errorInfo.message,'parameters differ');
            end
            checkpoint=d.checkpoint; checkpoint.version='TH_v14_20260930_joint_supply';
            oldFile=fullfile(folder,'old_checkpoint.mat'); save(oldFile,'checkpoint');
            o.initial.checkpointFile=oldFile; oldRejected=false;
            try
                TH_2D_segment_v15(o);
            catch errorInfo
                oldRejected=contains(errorInfo.message,'Incompatible checkpoint');
            end
            report=struct('timeDay',continuation.timeDay, ...
                'initialX',continuation.initialX,'checkpointX',d.checkpoint.X, ...
                'finalT',continuation.T,'referenceT',reference.T, ...
                'finalH',continuation.h,'referenceH',reference.h, ...
                'oldCheckpointRejected',oldRejected,'changedQuadratureRejected',changedRejected, ...
                'firstTime',first.timeDay);
        end
    end
end
