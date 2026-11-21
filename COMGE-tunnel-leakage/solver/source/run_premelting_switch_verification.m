function summary = run_premelting_switch_verification(folder)
%RUN_PREMELTING_SWITCH_VERIFICATION One full solved cooling/warming cycle.
% Three-domain closure is opt-in. This is NOT the historical 20-year case.
arguments
    folder (1,:) char
end
assert(~isfolder(folder),'PhaseSwitch:OutputExists','Preserve existing evidence; select a fresh folder.');
mkdir(folder); source=fullfile(folder,'source'); mkdir(source);
files=dir('*.m'); manifest=cell(numel(files),2);
for k=1:numel(files)
    copyfile(files(k).name,source);
    bytes=fileread(files(k).name); md=java.security.MessageDigest.getInstance('SHA-256');
    hash=typecast(md.digest(unicode2native(bytes,'UTF-8')),'uint8');
    manifest(k,:)={files(k).name,lower(reshape(dec2hex(hash,2)',1,[]))};
end
save(fullfile(folder,'source_manifest.mat'),'manifest');
p=joint_threephase_parameters(); f=premelting_switch_parameters(); c=joint_cohesive_parameters(); c.enabled=false;
m=build_cohesive_threephase_ring(make_thm_verification_mesh(16),3,p,c);
base=struct('cycles',1,'dtDays',0.5,'sealed',true,'verbose',true);
base.saveFile=fullfile(folder,'sealed_half_day.mat');
half=solve_premelting_switch_ring(m,p,f,base);
base.dtDays=0.25; base.saveFile=fullfile(folder,'sealed_quarter_day.mat');
quarter=solve_premelting_switch_ring(m,p,f,base);
base.dtDays=0.25; base.sealed=false; base.saveFile=fullfile(folder,'open_quarter_day.mat');
open=solve_premelting_switch_ring(m,p,f,base);
base.dtDays=0.125; base.saveFile=fullfile(folder,'open_eighth_day.mat');
eighth=solve_premelting_switch_ring(m,p,f,base);
testNames={'TestJointThermodynamics.m','TestNonisothermalInterface.m','TestThreePhaseInterface.m','TestCohesiveThreePhase.m','TestPremelting.m','TestPremeltingSwitch.m'};
results=runtests(testNames);
testSummary=struct('count',numel(results),'passed',sum([results.Passed]), ...
    'failed',sum([results.Failed]),'incomplete',sum([results.Incomplete]));
save(fullfile(folder,'matlab_tests.mat'),'results','testSummary');
summary=struct('scope','Component three-domain switching, solved heat + elastic Q4 ring; no soil or calibrated damage', ...
    'halfSteps',numel(half.time_s)-1,'quarterSteps',numel(quarter.time_s)-1,'openSteps',numel(open.time_s)-1,'openFineSteps',numel(eighth.time_s)-1, ...
    'testsPassed',all([results.Passed]),'testsCount',numel(results),'independentAuditPending',true);
save(fullfile(folder,'campaign_summary.mat'),'summary');
end
