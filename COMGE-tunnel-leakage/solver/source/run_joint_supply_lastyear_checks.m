function study=run_joint_supply_lastyear_checks(studyFolder,caseNames)
%RUN_JOINT_SUPPLY_LASTYEAR_CHECKS V14 final-year timestep sensitivity.
% study=run_joint_supply_lastyear_checks(STUDYFOLDER) examines base and
% supplySeasonal. An optional cell/string list also permits supply and
% transition. Each source must contain its own complete 19- or 20-year
% v14 history, including annual diagnostics. A twenty-year source supplies
% lastYear.initialX; a nineteen-year source supplies lastYear.T/h.
%
% The source mesh, material parameters and forcing are held fixed. One-day
% and half-day year-20 runs start from the SAME end-of-year-19 state. A
% nineteen-year source also gets a four-day comparison run. This checks
% final-year integration only. A different mesh must start a new complete
% history; no coarse-grid field is transferred to a fine grid here.
root=fileparts(mfilename('fullpath')); addpath(root);
if nargin<1||isempty(studyFolder)
    studyFolder=fullfile(root,'joint_supply_study_20260930');
end
if nargin<2||isempty(caseNames), caseNames={'base','supplySeasonal'}; end
studyFolder=char(studyFolder); caseNames=cellstr(string(caseNames));
assert(all(ismember(caseNames,{'base','transition','supply','supplySeasonal'})), ...
    'Choose base, transition, supply or supplySeasonal sources.');
assert(numel(unique(caseNames))==numel(caseNames),'Case names cannot repeat.');
outputFolder=fullfile(studyFolder,'year20_timestep_checks');
if ~isfolder(outputFolder), mkdir(outputFolder); end
definition=struct('solverVersion','TH_v14_20260930_joint_supply', ...
    'startDay',19*365,'durationDays',365,'referenceMaximumStep_day',.5, ...
    'caseNames',{caseNames},'outputFolder',outputFolder, ...
    'limitation',['Only final-year time integration is compared. The previous ' ...
        'nineteen-year state and spatial mesh are fixed for each case. This ' ...
        'does not prove full-history time convergence, spatial convergence, ' ...
        'physical calibration or periodic steady state. Fine meshes require ' ...
        'their own fresh full histories.']);
comparison=table(); caseDefinitions=cell(numel(caseNames),1);
allOutputs=cell(numel(caseNames),1);
for caseIndex=1:numel(caseNames)
    name=caseNames{caseIndex};
    sourceFile=fullfile(studyFolder,name,'TH_segment_spinup_results.mat');
    assert(isfile(sourceFile),'Missing complete v14 source: %s.',sourceFile);
    loaded=load(sourceFile,'result'); r=loaded.result;
    [seed,existing,caseDefinition]=sourceState(r,sourceFile,name);
    caseDefinitions{caseIndex}=caseDefinition;
    outFolder=fullfile(outputFolder,name);
    if ~isfolder(outFolder), mkdir(outFolder); end
    outputs=cell(3,1); elapsed=nan(3,1); dtMaxDays=[4;1;.5];
    labels={'dt4';'dt1';'dtHalf'};
    firstRun=1;
    if ~isempty(existing)
        outputs{1}=existing;
        dtMaxDays(1)=existing.parameters.dtMax/86400;
        labels{1}='existing_full_history'; firstRun=2;
        assertCompleteYear(existing);
        exportCase(existing,outFolder,labels{1});
    end
    for index=firstRun:3
        opt=seed.parameters; opt.mesh=seed.mesh;
        opt.initial=struct('mode','fields','T',seed.T,'h',seed.h,'startDay',19*365);
        opt.durationDays=365; opt.outputEveryDays=365/12;
        opt.dtMax=dtMaxDays(index)*86400;
        opt.dtInitial=min(.1*86400,opt.dtMax);
        opt.mechanics.enable=false; opt.plot=false; opt.verbose=false; opt.saveOutput=true;
        opt.outputPrefix=fullfile(outFolder,labels{index});
        fprintf('V14 year-20 check: %s/%s, same own year-19 state and mesh.\n',name,labels{index});
        started=tic; out=TH_2D_segment_v14(opt); elapsed(index)=toc(started);
        assertCompleteYear(out); outputs{index}=out;
        metadata=struct('caseName',name,'integrationLabel',labels{index}, ...
            'elapsed_s',elapsed(index),'caseDefinition',caseDefinition, ...
            'studyDefinition',definition); %#ok<NASGU>
        save([opt.outputPrefix '_results.mat'],'metadata','-append');
        exportCase(out,outFolder,labels{index});
    end
    Q=zeros(3,1); daily=cell(3,1); monthly=cell(3,1); rows=struct([]);
    for index=1:3
        out=outputs{index}; h=out.history;
        Q(index)=1000*out.leakage.totalVolume_m3_per_m;
        daily{index}=1000*out.leakage.dailyTotalVolume_m3_per_m(:);
        monthly{index}=1000*out.leakage.monthlyTotalVolume_m3_per_m(:);
        unrestricted=1000*sum((h(:,2)*86400).* ...
            sum(out.leakage.hydraulics.unrestrictedRate_m2_per_s,2));
        row=struct('caseName',string(name),'integrationLabel',string(labels{index}), ...
            'maximumStep_day',dtMaxDays(index),'Q20_L_per_m',Q(index), ...
            'sameTemperatureHydrostaticSeasonal_L_per_m',unrestricted, ...
            'sourceCompletedYears',r.completedYears, ...
            'meshRefinementLevel',caseDefinition.meshRefinementLevel, ...
            'thermalNodes',numel(out.T),'soilHeadNodes',numel(out.h), ...
            'acceptedSteps',out.solverWork.acceptedSteps, ...
            'rejectedSteps',out.solverWork.rejectedSteps,'elapsed_s',elapsed(index), ...
            'minimumStep_day',min(h(:,2)),'maximumScaledResidual',max(h(:,11)), ...
            'signedWaterDefect_m2',sum(h(:,8)),'sumAbsWaterDefect_m2',sum(abs(h(:,8))), ...
            'maxAbsWaterDefect_m2',max(abs(h(:,8))), ...
            'signedEnergyDefect_J_per_m',sum(h(:,9)), ...
            'sumAbsEnergyDefect_J_per_m',sum(abs(h(:,9))), ...
            'maxAbsEnergyDefect_J_per_m',max(abs(h(:,9))), ...
            'maximumNegativeTransEdges',max(out.leakage.hydraulics.negativeTransEdges));
        if index==1, rows=row; else, rows(index)=row; end %#ok<AGROW>
    end
    T=struct2table(rows);
    T.annualDifferenceRelativeToHalfDay_percent=100*abs(Q-Q(3))/max(abs(Q(3)),realmin);
    dailyDifference=zeros(3,1); monthlyDifference=zeros(3,1);
    for index=1:3
        dailyDifference(index)=sum(abs(daily{index}-daily{3}));
        monthlyDifference(index)=max(abs(monthly{index}-monthly{3}));
    end
    T.dailyL1DifferenceRelativeToHalfDay_L_per_m=dailyDifference;
    T.dailyL1RelativeToHalfDay_percent=100*dailyDifference/max(abs(Q(3)),realmin);
    T.maximumMonthlyDifferenceRelativeToHalfDay_L_per_m=monthlyDifference;
    writetable(T,fullfile(outFolder,'year20_timestep_comparison.csv'));
    if isempty(comparison), comparison=T; else, comparison=[comparison;T]; end %#ok<AGROW>
    allOutputs{caseIndex}=outputs;
end
study=struct('definition',definition,'caseDefinitions',{caseDefinitions}, ...
    'comparison',comparison,'outputs',{allOutputs});
writetable(comparison,fullfile(outputFolder,'year20_timestep_comparison.csv'));
save(fullfile(outputFolder,'year20_timestep_study.mat'),'study','-v7.3');
disp(comparison(:,{'caseName','integrationLabel','Q20_L_per_m', ...
    'annualDifferenceRelativeToHalfDay_percent','dailyL1RelativeToHalfDay_percent', ...
    'rejectedSteps','elapsed_s'}));
fprintf('%s\n',definition.limitation);
end

function [seed,existing,definition]=sourceState(r,sourceFile,name)
assert(isstruct(r)&&isfield(r,'completedYears')&&ismember(r.completedYears,[19 20]) ...
    &&isfield(r,'lastYear')&&strcmp(r.lastYear.version,'TH_v14_20260930_joint_supply'), ...
    'Source must be its own complete nineteen- or twenty-year v14 result.');
assert(isfield(r,'diagnostics')&&isfield(r,'columns'), ...
    'Source lacks its complete annual diagnostic history.');
yearColumn=find(strcmp(r.columns,'year'),1);
assert(~isempty(yearColumn)&&isequal(r.diagnostics(:,yearColumn),(1:r.completedYears)'), ...
    'Source annual diagnostics must contain every physical year from one.');
out=r.lastYear; nT=numel(out.T); existing=[];
assert(abs(out.timeDay-r.completedYears*365)<1e-7, ...
    'Source does not end at its stated completed year.');
if r.completedYears==20
    assert(numel(out.initialX)==nT+numel(out.h) ...
        &&abs(out.history(1,1)-out.history(1,2)-19*365)<1e-7, ...
        'The source lacks its actual end-of-year-19 initial state.');
    T=out.initialX(1:nT); h=out.initialX(nT+1:end); existing=out;
    initialSource='own complete year-20 lastYear.initialX (end of year 19)';
else
    assert(abs(sum(out.history(:,2))-365)<1e-6 ...
        &&abs(out.history(1,1)-out.history(1,2)-18*365)<1e-7, ...
        'The nineteen-year source has an incomplete final year.');
    T=out.T; h=out.h;
    initialSource='own complete year-19 lastYear.T/h';
end
assert(all(isfinite(T))&&all(isfinite(h)),'Year-19 fields contain nonfinite values.');
level=NaN;
if isfield(out.mesh,'jointRefinementReport'), level=out.mesh.jointRefinementReport.level; end
seed=struct('T',T,'h',h,'parameters',out.parameters,'mesh',out.mesh);
definition=struct('caseName',name,'sourceFile',sourceFile, ...
    'sourceCompletedYears',r.completedYears,'initialStateSource',initialSource, ...
    'meshRefinementLevel',level,'stateStartDay',19*365, ...
    'parameterPolicy','all physical parameters and forcing copied from this same source', ...
    'meshPolicy','same source mesh, no field transfer or mesh substitution');
end

function assertCompleteYear(out)
h=out.history; leak=out.leakage; nj=numel(out.parameters.leakage.jointAnglesRad);
assert(strcmp(out.version,'TH_v14_20260930_joint_supply') ...
    &&abs(out.timeDay-20*365)<1e-7&&abs(sum(h(:,2))-365)<1e-6 ...
    &&all(isfinite(h(:)))&&all(diff(h(:,1))>0), ...
    'The comparison did not finish a complete accepted twentieth year.');
assert(isequal(size(leak.jointFreeze.history),[size(h,1),4*nj]) ...
    &&isequal(size(leak.hydraulics.history),[size(h,1),3*nj+3]), ...
    'Thermal or hydraulic joint histories are incomplete.');
assert(all(isfinite(leak.hydraulics.history(:))), ...
    'Hydraulic diagnostic history contains nonfinite values.');
assert(out.parameters.liningK==0&&leak.bodyVolume_m3_per_m==0 ...
    &&abs(sum(leak.volumeByJoint_m3_per_m)-leak.totalVolume_m3_per_m)<1e-9, ...
    'The conservative joint-only volume partition failed.');
assert(numel(leak.dailyTotalVolume_m3_per_m)==365 ...
    &&numel(leak.monthlyTotalVolume_m3_per_m)==12, ...
    'Daily or monthly volume integration is incomplete.');
end

function exportCase(out,folder,label)
H=array2table(out.history,'VariableNames',out.historyColumns);
freeze=out.leakage.jointFreeze; hydraulic=out.leakage.hydraulics;
nj=numel(out.parameters.leakage.jointAnglesRad);
for joint=1:nj
    H.(sprintf('joint%d_factor',joint))=freeze.factor(:,joint);
    H.(sprintf('joint%d_drivingHead_m',joint))=hydraulic.drivingHead_m(:,joint);
    H.(sprintf('joint%d_referenceHead_m',joint))=hydraulic.referenceHead_m(:,joint);
    H.(sprintf('joint%d_unrestrictedRate_m2_per_s',joint))=hydraulic.unrestrictedRate_m2_per_s(:,joint);
end
H.waterTableDepth_m=hydraulic.waterTableDepth_m;
H.minimumTransmissivity_m_per_s=hydraulic.minTrans_m_per_s;
H.negativeTransmissivityEdges=hydraulic.negativeTransEdges;
writetable(H,fullfile(folder,[label '_accepted_history.csv']));
phaseDay=(1:365)'; dailyVolume_L_per_m=1000*out.leakage.dailyTotalVolume_m3_per_m(:);
writetable(table(phaseDay,dailyVolume_L_per_m),fullfile(folder,[label '_daily_volume.csv']));
phaseMonthIndex=(1:12)'; monthlyVolume_L_per_m=1000*out.leakage.monthlyTotalVolume_m3_per_m(:);
writetable(table(phaseMonthIndex,monthlyVolume_L_per_m),fullfile(folder,[label '_monthly_volume.csv']));
end
