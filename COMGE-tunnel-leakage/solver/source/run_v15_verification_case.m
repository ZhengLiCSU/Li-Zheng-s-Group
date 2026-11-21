function report=run_v15_verification_case(caseName)
%RUN_V15_VERIFICATION_CASE Independent first-year CVFE verification.
% Existing v14 histories and outputs are never reused or overwritten.
if nargin<1, caseName='annual'; end
root=fileparts(mfilename('fullpath')); addpath(root); maxNumCompThreads(4);
folder=fullfile(root,'v15_verification_20260930',caseName);
assert(~isfile(fullfile(folder,'verification.mat')),'Verification case already exists.');
if ~isfolder(folder), mkdir(folder); end
o=joint_supply_th_model(struct('durationDays',365,'dtMax',86400, ...
    'outputEveryDays',365/12,'plot',false,'saveOutput',true,'verbose',true, ...
    'outputPrefix',fullfile(folder,'cvfe'),'mechanics',struct('enable',false)));
switch caseName
    case 'annual'
        o.integration=struct('soilQuadratureOrder',5,'subdivisionLevel',0);
    case 'annualQuadrature2'
        o.integration=struct('soilQuadratureOrder',2,'subdivisionLevel',0);
    otherwise
        error('V15Verification:Case','Unknown case %s.',caseName);
end
started=tic; out=TH_2D_segment_v15(o); elapsed=toc(started);
report=struct('caseName',caseName,'elapsed_s',elapsed, ...
    'version',out.version,'durationDays',out.timeDay, ...
    'leak_L_per_m',1000*out.leakage.totalVolume_m3_per_m, ...
    'monthlyLeak_L_per_m',1000*out.leakage.monthlyTotalVolume_m3_per_m, ...
    'maximumWaterDefect_m2',max(abs(out.history(:,8))), ...
    'maximumEnergyDefect_Jm',max(abs(out.history(:,9))), ...
    'acceptedSteps',out.solverWork.acceptedSteps,'rejectedSteps',out.solverWork.rejectedSteps, ...
    'minimumStep_day',min(out.history(:,2)), ...
    'quadratureOrder',out.parameters.integration.soilQuadratureOrder);
assert(out.timeDay==365&&out.leakage.bodyVolume_m3_per_m==0);
assert(report.maximumWaterDefect_m2<1e-6&&report.maximumEnergyDefect_Jm<100);
assert(all(isfinite(out.T))&&all(isfinite(out.h)));
save(fullfile(folder,'verification.mat'),'report','-v7');
disp(report);
end
