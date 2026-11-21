function result=run_joint_supply_fine_mesh_case(caseName)
%RUN_JOINT_SUPPLY_FINE_MESH_CASE Complete-history active-joint mesh check.
% With no argument or 'supplySeasonal', repeats supplySeasonal from the
% initial profile on a separately generated level-2 mesh: Ksat=1e-10 m/s,
% joint-7 Cref=2e-10 m/s, and water-table depth 5 +/- 0.5 m with peak head
% at day 135. With 'base', repeats Ksat=1e-7 m/s and fixed 5 m water table.
% No coarse-mesh field is interpolated into this twenty-year calculation.
% Outputs are kept in supplySeasonalFine or baseFine, respectively, under
% joint_supply_study_20260930. Resume only the same case after interruption.
if nargin<1||isempty(caseName), caseName='supplySeasonal'; end
caseName=char(caseName);
assert(any(strcmp(caseName,{'base','supplySeasonal'})), ...
 'Select base or supplySeasonal for the complete-history mesh check.');
seasonal=strcmp(caseName,'supplySeasonal');
fineName=[caseName 'Fine'];
root=fileparts(mfilename('fullpath')); addpath(root);
folder=fullfile(root,'joint_supply_study_20260930',fineName);
if ~isfolder(folder), mkdir(folder); end
base=joint_freeze_th_model();
active=find(base.leakage.jointConductance>0);
assert(isequal(active,7),'Mesh verification requires the illustrative joint-7 defect.');
width=base.leakage.jointWidth_m;
if ~isscalar(width), width=width(active); end
meshFile=fullfile(folder,'mesh_active_joint7_refined_level2.mat');
reportFile=fullfile(folder,'mesh_refinement_report.mat');
if ~isfile(meshFile)
 S=load(fullfile(root,'example_mesh_repaired_interface.mat'));
 [Sfine,report]=refine_active_joint_mesh(S,base.leakage.jointAnglesRad(active),width,2);
 save(meshFile,'-struct','Sfine','-v7');
 save(reportFile,'report','-v7');
else
 assert(isfile(reportFile),'Fine mesh lacks its generation report.');
 d=load(reportFile,'report'); report=d.report;
end
Ksat=1e-7; amplitude=0; peakDay=NaN;
if seasonal, Ksat=1e-10; amplitude=.5; peakDay=135; end
model=struct('meshFile',meshFile,'Ksat',Ksat,'dtMax',4*86400, ...
 'initial',struct('mode','profile','startDay',0), ...
 'leakage',struct('jointConductance',base.leakage.jointConductance));
if seasonal
 day=(0:365)'; depth=5-.5*cos(2*pi*(day-135)/365); depth(end)=depth(1);
 model.forcing=struct('day',day,'waterTableDepth_m',depth);
end
caseInfo=struct('caseName',fineName, ...
 'purpose',['fresh full-history comparison against ' caseName ' level 1'], ...
 'Ksat_m_per_s',Ksat,'joint7_Cref_m_per_s',base.leakage.jointConductance(7), ...
 'waterTableMeanDepth_m',5,'waterTableAmplitude_m',amplitude, ...
 'peakExternalHeadPhaseDay',peakDay,'meshRefinementLevel',2, ...
 'meshFile',meshFile,'meshRefinementReport',report, ...
 'initialHistory','fresh profile, own full 20 years', ...
 'solverVersion','TH_v14_20260930_joint_supply','threads',4,'maximumStep_day',4);
save(fullfile(folder,'case_config.mat'),'caseInfo','model','-v7');
previousThreads=maxNumCompThreads(4);
cleanThreads=onCleanup(@() maxNumCompThreads(previousThreads)); %#ok<NASGU>
try
 result=run_TH_segment_year20_regime(struct('model',model,'outputFolder',folder, ...
  'resume',isfile(fullfile(folder,'spinup_progress.mat'))));
 assert(result.completedYears==20&& ...
  strcmp(result.lastYear.version,'TH_v14_20260930_joint_supply'), ...
  'Fine-mesh check requires a completed v14 twenty-year history.');
 result.studyCase=caseInfo;
 save(fullfile(folder,'TH_segment_spinup_results.mat'),'result','-v7');
 fid=fopen(fullfile(folder,'COMPLETE.txt'),'w');
 assert(fid>=0,'Cannot write fine-mesh completion marker.');
 cleanFile=onCleanup(@() fclose(fid)); %#ok<NASGU>
 fprintf(fid,'case=%s\nQ20_L_per_m=%.12g\nruntime_s=%.9g\n', ...
  fineName,result.year20AnnualLeakage_L_per_m,result.elapsedWallTime_s);
catch ME
 fid=fopen(fullfile(folder,'ERROR.txt'),'w');
 if fid>=0, fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid); end
 rethrow(ME);
end
end
