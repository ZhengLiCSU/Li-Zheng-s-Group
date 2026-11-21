function [result,caseInfo]=launch_joint_supply_study(caseName)
%LAUNCH_JOINT_SUPPLY_STUDY Independent complete-history supply-limit cases.
% launch_joint_supply_study('base'), ('transition'), ('supply'),
% ('supplySeasonal') or ('independentNoDrain') may run in separate MATLAB
% processes. With no argument, all five cases run sequentially.
% Every case starts from the same profile on its own level-1 refined mesh.
% No case borrows a field or checkpoint from the legacy v12 baseline.
%
% The seasonal water-table scenario has mean depth 5 m and amplitude 0.5 m.
% Its shallowest level (largest external water head) is on phase day 135,
% approximately mid-May in the repeated non-leap year. It is hypothetical.
% IndependentNoDrain solves the background TH history without tunnel water
% extraction, supplying an independent thermal-input leakage approximation.
if nargin<1 || isempty(caseName), caseName='all'; end
caseName=char(caseName);
root=fileparts(mfilename('fullpath')); addpath(root);
names={'base','transition','supply','supplySeasonal','independentNoDrain'};
if strcmp(caseName,'all')
 result=struct(); caseInfo=struct();
 for k=1:numel(names)
  [result.(names{k}),caseInfo.(names{k})]=launch_joint_supply_study(names{k});
 end
 return
end
assert(any(strcmp(caseName,names)), ...
 'Select base, transition, supply, supplySeasonal, independentNoDrain or all.');
studyFolder=fullfile(root,'joint_supply_study_20260930');
folder=fullfile(studyFolder,caseName);
if ~isfolder(folder), mkdir(folder); end
base=joint_freeze_th_model();
active=find(base.leakage.jointConductance>0);
assert(isequal(active,7),'The paired study requires the illustrative crown defect at joint 7.');
width=base.leakage.jointWidth_m;
if ~isscalar(width), width=width(active); end
meshFile=fullfile(folder,'mesh_active_joint7_refined_level1.mat');
meshReportFile=fullfile(folder,'mesh_refinement_report.mat');
if ~isfile(meshFile)
 S=load(fullfile(root,'example_mesh_repaired_interface.mat'));
 [Sfine,report]=refine_active_joint_mesh(S,base.leakage.jointAnglesRad(active),width,1);
 % The engine loads mesh MAT fields directly; preserve the full mesh struct.
 save(meshFile,'-struct','Sfine','-v7');
 save(meshReportFile,'report','-v7');
else
 assert(isfile(meshReportFile),'The refined mesh lacks its generation report.');
 d=load(meshReportFile,'report'); report=d.report;
end
Ksat=1e-7;
if strcmp(caseName,'transition'), Ksat=1e-9; end
if any(strcmp(caseName,{'supply','supplySeasonal','independentNoDrain'})), Ksat=1e-10; end
seasonal=any(strcmp(caseName,{'supplySeasonal','independentNoDrain'}));
model=struct('meshFile',meshFile,'Ksat',Ksat,'dtMax',4*86400, ...
 'initial',struct('mode','profile','startDay',0), ...
 'leakage',struct('jointConductance',base.leakage.jointConductance));
amplitude=0; peakHeadDay=NaN;
if seasonal
 amplitude=.5; peakHeadDay=135;
 day=(0:365)'; depth=5-amplitude*cos(2*pi*(day-peakHeadDay)/365);
 depth(end)=depth(1); % exact repeated endpoint for forcing validation
 model.forcing=struct('day',day,'waterTableDepth_m',depth);
end
if strcmp(caseName,'independentNoDrain'), model.tunnelDrain=false; end
caseInfo=struct('caseName',caseName,'purpose','uncalibrated parameterized mechanism study', ...
 'Ksat_m_per_s',Ksat,'joint7_Cref_m_per_s',base.leakage.jointConductance(7), ...
 'waterTableMeanDepth_m',5,'waterTableAmplitude_m',amplitude, ...
 'peakExternalHeadPhaseDay',peakHeadDay,'meshRefinementLevel',1, ...
 'meshFile',meshFile,'meshRefinementReport',report, ...
 'initialHistory','fresh profile, own full 20 years', ...
 'solverVersion','TH_v14_20260930_joint_supply', ...
 'noTunnelExtraction',strcmp(caseName,'independentNoDrain'), ...
 'threads',4,'maximumStep_day',4);
save(fullfile(folder,'case_config.mat'),'caseInfo','model','-v7');
previousThreads=maxNumCompThreads(4);
cleanThreads=onCleanup(@() maxNumCompThreads(previousThreads)); %#ok<NASGU>
user=struct('model',model,'outputFolder',folder, ...
 'resume',isfile(fullfile(folder,'spinup_progress.mat')));
try
 result=run_TH_segment_year20_regime(user);
 assert(result.completedYears==20,'The case did not finish its full 20-year history.');
 assert(strcmp(result.lastYear.version,'TH_v14_20260930_joint_supply'), ...
  'The case must be produced by the v14 solver.');
 result.studyCase=caseInfo;
 save(fullfile(folder,'TH_segment_spinup_results.mat'),'result','-v7');
 fid=fopen(fullfile(folder,'COMPLETE.txt'),'w');
 assert(fid>=0,'Cannot write case completion marker.');
 closeFile=onCleanup(@() fclose(fid)); %#ok<NASGU>
 fprintf(fid,'case=%s\nQ20_L_per_m=%.12g\nruntime_s=%.9g\n', ...
  caseName,result.year20AnnualLeakage_L_per_m,result.elapsedWallTime_s);
catch ME
 fid=fopen(fullfile(folder,'ERROR.txt'),'w');
 if fid>=0
  fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
 end
 rethrow(ME);
end
end
