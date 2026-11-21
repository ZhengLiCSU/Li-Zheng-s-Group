function report=verify_thm20_time_step(source)
%VERIFY_THM20_TIME_STEP Target-year THM time-step check from identical Y19.
% REPORT=verify_thm20_time_step(FOLDER) reads a COMPLETED twenty-year
% run_TH_segment_year20_thm folder. Its returned RESULT struct is also an
% accepted input; RESULT.config.outputFolder must retain the core files.
% Default folder: thm20_20260930_analytic_p1 beside this function.
%
% The original year-20 maximum step must be one day. The independent fine
% calculation starts from the SAME complete year-19 THM checkpoint, with a
% half-day maximum step. Solver tolerances, all physical laws/parameters,
% fifth-order soil phase integration, mesh and installed U/T are preserved.
% This is TARGET-YEAR temporal convergence conditional on the original Y19
% history; it is NOT a comparison of two independently refined 20-year runs.
% No parameter changes or relaxed guards are applied if the check fails.
%
% Primary criterion: annual leakage difference / fine annual leakage <=1%.
% Monthly differences are also required to be <=1% of FINE ANNUAL leakage;
% ordinary monthly relative differences are reported without failing tiny
% frozen-month volumes. T/h/U/strain errors are diagnostics, without an
% arbitrary 1% field criterion. Original solver-derived budget bounds and
% actual model domain bounds are independently checked for both runs.
% A compatible complete fine calculation is reused, never an incomplete
% checkpoint or a fine result associated with a different starting state.

root=fileparts(mfilename('fullpath')); addpath(root);
if nargin<1||isempty(source), source=fullfile(root,'thm20_20260930_analytic_p1'); end
[result,folder,sourceFile]=readSource(source);
assert(isfield(result,'completedYears')&&result.completedYears==20 ...
    &&isfield(result,'engineVersion')&&strcmp(result.engineVersion,'THM_v16_20260930_analytic_p1'), ...
    'THM20Time:IncompleteSource','A completed twenty-year analytic-P1 THM result is required.');
yearColumn=find(strcmp(result.annualColumns,'year'),1);
assert(~isempty(yearColumn)&&size(result.annual,1)==20 ...
    &&isequal(result.annual(:,yearColumn),(1:20)')&&numel(result.coreCheckpoints)==20, ...
    'THM20Time:History','Source must register all twenty consecutive completed THM years.');
reference=result.lastYear; checkedYear(reference);
p=reference.parameters;
assert(abs(p.dtMax-86400)<1e-9&&p.integration.soilQuadratureOrder==5, ...
    'THM20Time:Discretization','This check requires the original one-day maximum and fifth-order phase quadrature.');
cp19File=checkpointPath(result.coreCheckpoints{19},folder,19);
cp20File=checkpointPath(result.coreCheckpoints{20},folder,20);
cp19=readCheckpoint(cp19File,19*365); cp20=readCheckpoint(cp20File,20*365);
installed=result.installedMechanicalReference;
assert(isequal(cp19.X,reference.initialX)&&isequal(cp20.X,reference.finalX) ...
    &&isequal(cp19.mechanicalReference,installed)&&isequal(cp20.mechanicalReference,installed), ...
    'THM20Time:Reference','Source checkpoints must match full initial/final X and the unchanged installed reference.');
assert(isequal(cp19.mechanicalReference.U,reference.mechanics.referenceU) ...
    &&isequal(cp19.mechanicalReference.T,reference.mechanics.referenceTemperature) ...
    &&isequaln(cp19.signature,physicalSignature(p)) ...
    &&isequal(cp19.pAll,reference.mesh.pAll)&&isequal(cp19.tSoil,reference.mesh.tSoil), ...
    'THM20Time:SourcePhysics','Source checkpoint has a different mesh, physical law or installed U/T.');

outputFolder=fullfile(folder,'year20_time_step_verification');
if ~isfolder(outputFolder), mkdir(outputFolder); end
prefix=fullfile(outputFolder,'year20_dtHalf');
fineFile=[prefix '_results.mat'];
opt=p; opt.mesh=reference.mesh;
opt.initial=struct('mode','checkpoint','checkpointFile',cp19File);
opt.dtMax=.5*86400; opt.dtInitial=min(p.dtInitial,opt.dtMax);
opt.durationDays=365; opt.outputPrefix=prefix;
opt.plot=false; opt.verbose=false; opt.saveOutput=true;
assert(opt.dtMin<=opt.dtMax,'THM20Time:Discretization','Original minimum step exceeds the fine maximum.');
reuseKey=struct('schemaVersion',1,'engineVersion',reference.version, ...
    'year19X',cp19.X,'installedReference',installed,'mesh',reference.mesh, ...
    'referenceInitialX',reference.initialX,'referenceFinalX',reference.finalX, ...
    'referenceParameters',reference.parameters, ...
    'referenceAnnualLeakage',reference.leakage.totalVolume_m3_per_m, ...
    'referenceMonthlyLeakage',reference.leakage.monthlyTotalVolume_m3_per_m, ...
    'fineParameters',opt);
reused=false; recoveredMetadata=false; writeMetadata=false;
if isfile(fineFile)
    saved=load(fineFile,'out','fineMetadata');
    assert(isfield(saved,'out'),'THM20Time:CacheMismatch', ...
        'Existing fine result lacks its core output; it has been preserved.');
    fine=saved.out; reused=true;
    if isfield(saved,'fineMetadata')
        assert(isfield(saved.fineMetadata,'reuseKey') ...
            &&sameReuseKey(saved.fineMetadata.reuseKey,reuseKey), ...
            'THM20Time:CacheMismatch','Existing fine result belongs to a different source/state/settings; it has been preserved.');
        elapsed=saved.fineMetadata.elapsedWallTime_s;
    else
        % Core saves the complete result before this wrapper appends metadata.
        % An interruption in that interval must not strand a valid fine run.
        % Recovery is allowed ONLY after the independent full checks below.
        elapsed=NaN; recoveredMetadata=true; writeMetadata=true;
    end
else
    fprintf('Recomputing THM year 20 from complete Y19 checkpoint; dtMax=0.5 day, same physics/tolerances/quadrature.\n');
    started=tic; fine=TH_2D_segment_v16(opt); elapsed=toc(started);
    writeMetadata=true;
end
checkedYear(fine);
fineCP=readCheckpoint([prefix '_checkpoint.mat'],20*365);
assert(isequal(fine.initialX,cp19.X)&&isequal(fineCP.X,fine.finalX) ...
    &&isequal(fineCP.mechanicalReference,installed) ...
    &&isequaln(fineCP.signature,physicalSignature(p)) ...
    &&isequaln(fineCP.parameters,fine.parameters) ...
    &&isequal(fine.mechanics.referenceU,installed.U) ...
    &&isequal(fine.mechanics.referenceTemperature,installed.T), ...
    'THM20Time:FineReference','Fine run changed the inherited X or installed mechanical reference.');
assert(isequaln(physicalSignature(fine.parameters),physicalSignature(p)) ...
    &&isequaln(fine.parameters.solver,p.solver) ...
    &&isequaln(fine.mesh,reference.mesh)&&fine.parameters.dtMax==.5*86400 ...
    &&fine.parameters.dtInitial==opt.dtInitial&&fine.parameters.durationDays==365 ...
    &&fine.parameters.dtMin==p.dtMin ...
    &&fine.parameters.outputEveryDays==p.outputEveryDays ...
    &&strcmp(fine.parameters.initial.mode,'checkpoint') ...
    &&(~isfield(fine.parameters.initial,'resumeHistory')||~fine.parameters.initial.resumeHistory), ...
    'THM20Time:FinePhysics','Fine run changed physical parameters, solver settings, integration, mesh or reporting interval.');
if writeMetadata
    fineMetadata=struct('reuseKey',reuseKey,'elapsedWallTime_s',elapsed, ...
        'recoveredFromIndependentlyCheckedCompleteOutput',recoveredMetadata);
    save(fineFile,'fineMetadata','-append');
end
if reused, fprintf('Reusing independently checked complete year-20 half-day THM result.\n'); end

q0=reference.leakage.totalVolume_m3_per_m;
q1=fine.leakage.totalVolume_m3_per_m;
annualRelative=relativeDifference(q0,q1);
m0=reference.leakage.monthlyTotalVolume_m3_per_m(:);
m1=fine.leakage.monthlyTotalVolume_m3_per_m(:);
monthlyRelative=relativeDifference(m0,m1);
monthlyAnnualRelative=ratioTo(abs(m0-m1),q1);
monthly=table((1:12)',1000*m0,1000*m1,1000*abs(m0-m1), ...
    100*monthlyRelative,100*monthlyAnnualRelative, ...
    'VariableNames',{'month','reference_L_per_m','fine_L_per_m', ...
    'absoluteDifference_L_per_m','differenceRelativeToFineMonth_percent', ...
    'differenceRelativeToFineAnnual_percent'});
diagnostics=fieldDifferences(reference,fine);
[referenceChecks,referenceChecksPassed]=physicalChecks(reference,'reference');
[fineChecks,fineChecksPassed]=physicalChecks(fine,'fine');
checks=[referenceChecks;fineChecks];
primaryPassed=annualRelative<=.01&&all(monthlyAnnualRelative<=.01);
passed=primaryPassed&&referenceChecksPassed&&fineChecksPassed;
if passed, status='passed_target_year_leakage_time_step_check';
else, status='not_passed_results_retained'; end
scope=['Year-20 time-step convergence from the SAME complete year-19 THM state, ' ...
    'with the original installed mechanical reference and preceding nineteen-year history. ' ...
    'This is not convergence of two complete twenty-year histories, spatial convergence, ' ...
    'physical calibration, or validation of omitted segregated-ice/density-expansion mechanisms. ' ...
    'Passed refers to annual/monthly leakage criteria and the explicit numerical/domain checks; ' ...
    'field and strain errors remain diagnostics.'];
report=struct('schemaVersion',1,'status',status,'passed',passed, ...
    'primaryLeakagePassed',primaryPassed,'numericalAndDomainPassed',referenceChecksPassed&&fineChecksPassed, ...
    'sourceCompletedYears',20,'sourceFile',sourceFile,'sourceFolder',folder, ...
    'year19Checkpoint',cp19File,'fineResultFile',fineFile,'fineCheckpoint',[prefix '_checkpoint.mat'], ...
    'referenceDtMax_day',p.dtMax/86400,'fineDtMax_day',.5, ...
    'referenceAnnualLeakage_L_per_m',1000*q0,'fineAnnualLeakage_L_per_m',1000*q1, ...
    'annualAbsoluteDifference_L_per_m',1000*abs(q0-q1), ...
    'annualRelativeDifference_percent',100*annualRelative, ...
    'maximumMonthlyDifferenceRelativeToFineAnnual_percent',100*max(monthlyAnnualRelative), ...
    'leakageAcceptanceLimit_percent',1,'monthly',monthly,'fieldDiagnostics',diagnostics, ...
    'numericalAndDomainChecks',checks,'fineElapsedWallTime_s',elapsed, ...
    'compatibleFineResultReused',reused,'fineMetadataRecovered',recoveredMetadata,'scope',scope);
reportFile=fullfile(outputFolder,'year20_time_step_report.mat');
save(reportFile,'report','-v7');
writetable(monthly,fullfile(outputFolder,'year20_monthly_time_step_comparison.csv'));
writetable(diagnostics,fullfile(outputFolder,'year20_final_field_time_step_comparison.csv'));
writetable(checks,fullfile(outputFolder,'year20_numerical_domain_checks.csv'));
annual=table(1000*q0,1000*q1,1000*abs(q0-q1),100*annualRelative, ...
    100*max(monthlyAnnualRelative),primaryPassed,referenceChecksPassed&&fineChecksPassed,passed, ...
    'VariableNames',{'reference_L_per_m','fine_L_per_m','absoluteDifference_L_per_m', ...
    'annualRelativeDifference_percent','maxMonthlyDifferenceRelativeToFineAnnual_percent', ...
    'leakageCriterionPassed','numericalAndDomainPassed','passed'});
writetable(annual,fullfile(outputFolder,'year20_annual_time_step_comparison.csv'));
disp(annual); fprintf('%s\n',scope);
if ~passed
    warning('THM20Time:NotPassed','Target-year time-step check did not pass. Fine results and all diagnostics have been preserved without parameter changes.');
end
end

function [result,folder,file]=readSource(source)
file='';
if isstruct(source)&&isscalar(source)
    result=source;
    assert(isfield(result,'config')&&isfield(result.config,'outputFolder'), ...
        'THM20Time:Source','Result must retain config.outputFolder.');
    folder=canonicalPath(result.config.outputFolder);
else
    assert((ischar(source)&&isrow(source))||(isstring(source)&&isscalar(source)), ...
        'THM20Time:Source','Supply a completed THM folder or annual RESULT struct.');
    folder=canonicalPath(source);
    % The year-commit inventory is authoritative if an interruption occurs
    % after committing Y20 but before publishing the aggregate result file.
    file=fullfile(folder,'thm_progress.mat');
    if ~isfile(file), file=fullfile(folder,'TH_segment_year20_thm_results.mat'); end
    assert(isfile(file),'THM20Time:Source','Missing completed annual result: %s.',file);
    saved=load(file);
    if isfield(saved,'progress')
        result=saved.progress;
    elseif isfield(saved,'result')
        result=saved.result;
    else
        error('THM20Time:Source','Annual MAT must contain progress or result.');
    end
end
assert(isfolder(folder),'THM20Time:Source','Source folder does not exist.');
end

function file=canonicalPath(file)
file=char(file); object=java.io.File(file);
% Java's startup directory need not follow a subsequent MATLAB cd.
if ~object.isAbsolute(), file=fullfile(pwd,file); object=java.io.File(file); end
file=char(object.getCanonicalPath());
end

function file=checkpointPath(registered,folder,year)
[~,name,extension]=fileparts(char(registered));
assert(strcmp([name extension],sprintf('year_%02d_checkpoint.mat',year)), ...
    'THM20Time:Checkpoint','Registered checkpoint filename does not match its complete year.');
local=fullfile(folder,[name extension]);
if isfile(local), file=canonicalPath(local); else, file=canonicalPath(registered); end
end

function yes=sameReuseKey(a,b)
% Case relocation changes file locations, not starting X, solver settings,
% actual mesh, reference, or physical parameters. Ignore ONLY file-location
% fields and unused checkpoint-initial settings; full initial X is in key.
if ~isstruct(a)||~isscalar(a)||~isequal(sort(fieldnames(a)),sort(fieldnames(b)))
    yes=false; return
end
for name={'referenceParameters','fineParameters'}
    if ~isstruct(a.(name{1}))||~isscalar(a.(name{1})), yes=false; return; end
    remove={'meshFile','outputPrefix','initial'};
    a.(name{1})=rmfield(a.(name{1}),intersect(remove,fieldnames(a.(name{1}))));
    b.(name{1})=rmfield(b.(name{1}),intersect(remove,fieldnames(b.(name{1}))));
end
yes=isequaln(a,b);
end

function cp=readCheckpoint(file,day)
assert(isfile(file),'THM20Time:Checkpoint','Missing complete THM checkpoint: %s.',file);
saved=load(file,'checkpoint'); cp=saved.checkpoint;
assert(strcmp(cp.version,'THM_v16_20260930_analytic_p1')&&isfield(cp,'mechanicalReference') ...
    &&all(isfinite(cp.X(:)))&&abs(cp.time/86400-day)<1e-7, ...
    'THM20Time:Checkpoint','Checkpoint has incompatible physics, missing installed reference, or incorrect time.');
end

function checkedYear(out)
assert(strcmp(out.version,'THM_v16_20260930_analytic_p1')&&out.parameters.thm.enabled ...
    &&abs(out.timeDay-20*365)<1e-7&&isfield(out,'mechanics'), ...
    'THM20Time:Year','A complete year-20 two-way THM output is required.');
h=out.history;
assert(~isempty(h)&&all(isfinite(h(:)))&&all(diff(h(:,1))>0)&&all(h(:,2)>0) ...
    &&abs(sum(h(:,2))-365)<1e-6&&abs(h(1,1)-h(1,2)-19*365)<1e-7 ...
    &&all(abs(diff([19*365;h(:,1)])-h(:,2))<1e-6) ...
    &&size(out.mechanics.history,1)==size(h,1) ...
    &&size(out.mechanics.strainHistory,1)==size(h,1), ...
    'THM20Time:Year','The coupled accepted-step year-20 inventory is incomplete.');
assert(numel(out.leakage.monthlyTotalVolume_m3_per_m)==12 ...
    &&numel(out.leakage.dailyTotalVolume_m3_per_m)==365 ...
    &&abs(sum(out.leakage.monthlyTotalVolume_m3_per_m)-out.leakage.totalVolume_m3_per_m)<1e-8 ...
    &&abs(sum(out.leakage.dailyTotalVolume_m3_per_m)-out.leakage.totalVolume_m3_per_m)<1e-8 ...
    &&all(isfinite(out.finalX(:))), ...
    'THM20Time:Year','Target-year monthly/daily leakage or final coupled state is invalid.');
end

function p=physicalSignature(p)
p=thm_constitutive_signature(p);
end

function r=relativeDifference(a,b)
r=ratioTo(abs(a-b),abs(b));
end

function r=ratioTo(numerator,denominator)
if isscalar(denominator), denominator=repmat(denominator,size(numerator)); end
r=zeros(size(numerator)); positive=denominator>0;
r(positive)=numerator(positive)./denominator(positive);
r(~positive&numerator>0)=Inf;
end

function tab=fieldDifferences(a,b)
nT=numel(b.T); nh=numel(b.h); U0=a.finalX(nT+nh+1:end); U1=b.finalX(nT+nh+1:end);
refU=b.mechanics.referenceU; refT=b.mechanics.referenceTemperature;
base={a.T,a.T-refT,a.h,U0,U0-refU,a.mechanics.soilIncrementU_m, ...
    a.mechanics.ringU_m,a.mechanics.soilElementVolumetricStrain,a.mechanics.currentSoilPorosity};
fine={b.T,b.T-refT,b.h,U1,U1-refU,b.mechanics.soilIncrementU_m, ...
    b.mechanics.ringU_m,b.mechanics.soilElementVolumetricStrain,b.mechanics.currentSoilPorosity};
names={'finalTemperature_C','finalTemperatureIncrementFromInstalled_C','finalHead_m', ...
    'finalFullDisplacement_m','finalDisplacementIncrementFromInstalled_m', ...
    'finalSoilDisplacementIncrement_m','finalRingDisplacement_m', ...
    'finalSoilVolumetricStrain','finalSoilCurrentPorosity'};
absInf=zeros(numel(names),1); relL2=absInf;
for j=1:numel(names)
    difference=base{j}(:)-fine{j}(:);
    absInf(j)=norm(difference,Inf); relL2(j)=ratioTo(norm(difference),norm(fine{j}(:)));
end
% Values of max strain/contact quantities are compared over the complete
% target-year history, not at interpolated times or just final snapshots.
cols=[2 3 5 6 7]; sh0=a.mechanics.strainHistory; sh1=b.mechanics.strainHistory;
historyNames={'yearMaxIncrementalVolumetricStrain','yearMaxIncrementalPrincipalStrain', ...
    'yearMaxIncrementalGradientNorm','yearMaxAbsolutePrincipalStrain','yearMaxAbsoluteGradientNorm'};
for j=1:numel(cols)
    v0=max(sh0(:,cols(j))); v1=max(sh1(:,cols(j)));
    names{end+1}=historyNames{j}; absInf(end+1,1)=abs(v0-v1); %#ok<AGROW>
    relL2(end+1,1)=relativeDifference(v0,v1); %#ok<AGROW>
end
tab=table(string(names(:)),absInf,100*relL2, ...
    'VariableNames',{'quantity','absoluteInfinityDifference','relativeL2OrScalarDifference_percent'});
end

function [tab,passed]=physicalChecks(out,label)
p=out.parameters; h=out.history; s=out.mechanics.strainHistory;
i=out.mechanics.installedSoilStrainDiagnostics;
[vw,ve]=rowVolumes(out.mesh);
tol=p.solver.tol; waterBound=tol*p.solver.massScale*sum(vw);
energyBound=tol*p.solver.energyScale*sum(ve);
names={'maximumScaledResidual','maxLocalWaterDefect_m2','maxLocalEnergyDefect_J_per_m', ...
    'maxGlobalStepWaterDefect_m2','maxGlobalStepEnergyDefect_J_per_m', ...
    'annualAbsoluteWaterDefectSum_m2','annualAbsoluteEnergyDefectSum_J_per_m', ...
    'maxIncrementalVolumetricStrain','maxIncrementalPrincipalStrain','maxIncrementalGradientNorm', ...
    'maxAbsolutePrincipalStrain','maxAbsoluteGradientNorm','maxInstalledPrincipalStrain', ...
    'maxInstalledGradientNorm','maxOmittedNewSlitVolumeFraction','maxMechanicalForceResidual_N_per_m'};
values=[max(h(:,11));out.numerics.maximumAcceptedLocalWaterDefect_m2; ...
    out.numerics.maximumAcceptedLocalEnergyDefect_Jm;max(abs(h(:,8)));max(abs(h(:,9))); ...
    sum(abs(h(:,8)));sum(abs(h(:,9)));max(s(:,2));max(s(:,3));max(s(:,5));max(s(:,6));max(s(:,7)); ...
    i.maxAbsPrincipalStrain;i.maxGradientNorm;max(out.mechanics.history(:,6));max(out.mechanics.history(:,5))];
limits=[tol;tol*p.solver.massScale*max(vw);tol*p.solver.energyScale*max(ve); ...
    waterBound;energyBound;size(h,1)*waterBound;size(h,1)*energyBound; ...
    p.thm.maximumVolumetricStrain;p.thm.maximumPrincipalStrain;p.thm.maximumDisplacementGradient; ...
    p.thm.maximumPrincipalStrain;p.thm.maximumDisplacementGradient;p.thm.maximumPrincipalStrain; ...
    p.thm.maximumDisplacementGradient;p.thm.maxOmittedSlitVolumeFraction;tol/p.thm.forceResidualScale];
% Match validate_thm_20year's data-dependent floating-point allowance.
% The stored global inventory can be many orders larger than its defect;
% subtraction/summation roundoff is distinct from relaxing the Newton tol.
waterRound=128*eps(max(1,max(abs(h(:,4)))));
energyRound=128*eps(max(1,max(abs(h(:,5)))));
roundoff=zeros(size(limits));
roundoff(1)=tol*1e-7;
roundoff([2 4])=waterRound; roundoff([3 5])=energyRound;
roundoff(6)=size(h,1)*waterRound; roundoff(7)=size(h,1)*energyRound;
roundoff(16)=limits(16)*1e-6;
acceptance=limits+roundoff;
rowPassed=isfinite(values)&values>=0&values<=acceptance;
% These volume/energy bounds follow the original scaled residual equations;
% global sums are bounded by the sum of positive nodal row scales, without
% cancellation. No independently invented budget tolerance is used here.
comparison=repmat("<=",numel(names),1);
phi=out.mechanics.currentSoilPorosity;
names=[names,{'minimumFinalCurrentPorosity','maximumFinalCurrentPorosity'}];
values=[values;min(phi);max(phi)]; limits=[limits;0;1];
roundoff=[roundoff;0;0]; acceptance=[acceptance;0;1];
rowPassed=[rowPassed;isfinite(min(phi))&&min(phi)>0;isfinite(max(phi))&&max(phi)<1];
comparison=[comparison;">";"<"];
tab=table(repmat(string(label),numel(names),1),string(names(:)),values,comparison,limits,roundoff,acceptance,rowPassed, ...
    'VariableNames',{'run','check','observedMaximumOrSum','comparison', ...
    'modelOrSolverDerivedBound','floatingPointAllowance','actualAcceptanceBound','passed'});
passed=all(rowPassed)&&all(s(:,2:end)>=0,'all') ...
    &&all(out.mechanics.history(:,5:6)>=0,'all');
end

function [water,energy]=rowVolumes(S)
% Exact same positive T3 and 2x2 Q4 row volumes used by v16 residual scales.
N=size(S.pAll,1); water=zeros(N,1); concrete=water; t=S.tSoil;
a=S.pAll(t(:,1),:); b=S.pAll(t(:,2),:); c=S.pAll(t(:,3),:);
area=abs((b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(c(:,1)-a(:,1)).*(b(:,2)-a(:,2)))/2;
water=accumarray(t(:),repmat(area/3,3,1),[N 1]);
for name={'qTunnelGlobal','qInvertGlobal'}
    if ~isfield(S,name{1}), continue; end
    q=S.(name{1});
    for e=1:size(q,1)
        ids=q(e,:); x=S.pAll(ids,:); v=zeros(4,1);
        for xi=[-1 1]/sqrt(3)
            for eta=[-1 1]/sqrt(3)
                shape=.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
                d=.25*[-(1-eta),1-eta,1+eta,-(1+eta);-(1-xi),-(1+xi),1+xi,1-xi];
                v=v+shape'*abs(det(d*x));
            end
        end
        concrete(ids)=concrete(ids)+v;
    end
end
if isfield(S,'tInvertTriGlobal')&&~isempty(S.tInvertTriGlobal)
    t=S.tInvertTriGlobal; a=S.pAll(t(:,1),:); b=S.pAll(t(:,2),:); c=S.pAll(t(:,3),:);
    area=abs((b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(c(:,1)-a(:,1)).*(b(:,2)-a(:,2)))/2;
    concrete=concrete+accumarray(t(:),repmat(area/3,3,1),[N 1]);
end
energy=water+concrete; water=water(unique(S.tSoil(:)));
end
