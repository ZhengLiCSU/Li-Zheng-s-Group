function [report,figs]=validate_thm_20year(input,options)
%VALIDATE_THM_20YEAR Audit a complete, consecutive 20-year coupled history.
% [R,F]=validate_thm_20year(annualResultOrOutputFolder,OPTIONS) reads the
% annual result/progress MAT, EVERY year_XX_annual_results.mat, and EVERY
% registered year-end checkpoint. Tables and numeric matrices are accepted.
% Missing files/fields, nonfinite raw data and incomplete histories produce
% R.passed=false with explicit checks; they never become a successful audit.
% Supported annual schemas: 11/THM_v16_20260930_analytic_p1 and historical
% 10/THM_v16_20260930_segment_contact. A history cannot mix engine labels.
%
% OPTIONS: makePlots=true, saveReport=true, outputFolder='', visible='off',
% timeTolerance_day=1e-6, volumeAbsTolerance_m2=1e-8,
% relativeTolerance=1e-10, rateAbsTolerance_m2_s=1e-12.
% Plots are created only after a COMPLETE 20-year numerical audit passes.
% saveReport stores MAT/CSV/text and, if requested, scientific PNG/PDF plots
% under outputFolder (default: <case>/validation_20year). Existing raw case
% files are never modified. F is a struct of figure handles.
%
% r=validate_thm_20year('absolute/path/to/case');
% [r,f]=validate_thm_20year(result,struct('saveReport',true));
% Numerical consistency is not field/experimental validation. Tolerances
% below derive from the accepted nonlinear residual and nodal CVFE volumes;
% they do not establish physical accuracy or discretization convergence.
% The same-state simple flow uses the recorded TOTAL Ceff (capillary plus
% open slit), exact Robin quadrature and FAR head after the entry condition,
% not Cref times a freeze factor or positive part of an averaged head.

if nargin<2, options=struct(); end
o=struct('makePlots',true,'saveReport',true,'outputFolder','', ...
    'visible','off','timeTolerance_day',1e-6,'volumeAbsTolerance_m2',1e-8, ...
    'relativeTolerance',1e-10,'rateAbsTolerance_m2_s',1e-12);
assert(isstruct(options)&&isscalar(options),'THM20Audit:Options','Options must be a scalar struct.');
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'THM20Audit:Options','Unknown option: %s.',key{1});
    o.(key{1})=options.(key{1});
end
assert(isSwitch(o.makePlots)&&isSwitch(o.saveReport),'THM20Audit:Options','Plot/save switches must be logical scalars.');
assert(any(strcmp(o.visible,{'on','off'})),'THM20Audit:Options','visible must be on or off.');
for key={'timeTolerance_day','volumeAbsTolerance_m2','relativeTolerance','rateAbsTolerance_m2_s'}
    x=o.(key{1}); assert(isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0, ...
        'THM20Audit:Options','Tolerances must be positive finite scalars.');
end
figs=struct(); checks=struct('year',{},'name',{},'passed',{},'detail',{});
report=struct('validatorVersion','THM20Audit_20261002_cycle_memory','passed',false, ...
    'completedYears',NaN,'requiredYears',20,'schemaVersion',NaN,'engineVersion','', ...
    'inputFolder','','inputFile','','configuration',struct(),'parameters',struct(), ...
    'parameterSummary',struct(),'yearSummary',table(),'checks',table(), ...
    'issues',{{}},'warnings',{{}},'year20',struct(),'options',o, ...
    'sourceHash',sourceHashes(),'sourceAtStartHash',struct([]),'sources',researchSources(), ...
    'hydraulicAudit',struct([]),'numericalValidationOnly',true,'fieldValidation',false, ...
    'limitations',{limitations()},'plotSuppressionReason','', ...
    'createdAtLocal',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
folder=''; data=struct(); lastOut=[]; metrics=repmat(metricTemplate(),0,1); hydraulicAudits=struct([]);
try
    [data,folder,inputFile]=readInput(input);
    report.inputFolder=folder; report.inputFile=inputFile;
    if isfolder(fullfile(folder,'source_at_start'))
        report.sourceAtStartHash=sourceHashes(fullfile(folder,'source_at_start'));
    end
    need(data,{'schemaVersion','engineVersion','completedYears','config', ...
        'coreCheckpoints','installedMechanicalReference','lastYear'},'annual result/progress');
    report.schemaVersion=data.schemaVersion; report.engineVersion=char(data.engineVersion);
    report.configuration=data.config;
    check(0,'supported_schema_engine', ...
        (isequal(data.schemaVersion,11)&&strcmp(report.engineVersion,'THM_v16_20260930_analytic_p1')) ...
        ||(isequal(data.schemaVersion,10)&&strcmp(report.engineVersion,'THM_v16_20260930_segment_contact')), ...
        sprintf('schema=%s; engine=%s; no cross-version continuation accepted.',mat2str(data.schemaVersion),report.engineVersion));
    n=data.completedYears;
    demand(isnumeric(n)&&isreal(n)&&isscalar(n)&&isfinite(n)&&n>=0&&n<=20&&n==fix(n), ...
        'completedYears must be an integer from zero to twenty.');
    report.completedYears=n;
    check(0,'complete_twenty_years',n==20,sprintf('%d/20 complete years are registered.',n));
    check(0,'dedicated_case_folder',~isempty(folder)&&isfolder(folder), ...
        'All year archives and checkpoints are required, even when an annual result struct is supplied.');
    demand(iscell(data.coreCheckpoints),'coreCheckpoints must be a cell array.');
    check(0,'checkpoint_registration_count',numel(data.coreCheckpoints)==n, ...
        sprintf('%d checkpoints registered for %d years.',numel(data.coreCheckpoints),n));
    if n>0
        need(data.installedMechanicalReference,{'U','T'},'installed reference');
        demand(finite(data.installedMechanicalReference.U)&&finite(data.installedMechanicalReference.T), ...
            'Installed U/T reference must be finite numeric data.');
    end
    expected=columnSchema(); names=fieldnames(expected);
    ledger=struct();
    for k=1:numel(names)
        name=names{k}; meta=expected.(name).metadata;
        need(data,{name},'aggregate annual result/progress');
        if isfield(data,meta), columns=data.(meta); else, columns={}; end
        ledger.(name)=asMatrix(data.(name),columns,expected.(name).columns,name);
    end
    check(0,'annual_year_inventory',isequal(ledger.annual(:,1),(1:n)'), ...
        'Annual rows must contain every year once, in order.');
    check(0,'aggregate_row_counts',size(ledger.monthly,1)==12*n ...
        &&size(ledger.byJoint,1)==8*n&&size(ledger.daily,1)==365*n, ...
        sprintf('Annual/month/joint/day rows: %d/%d/%d/%d.',size(ledger.annual,1), ...
        size(ledger.monthly,1),size(ledger.byJoint,1),size(ledger.daily,1)));
    prior=[]; firstSignature=[]; firstGeometry=[];
    for year=1:n
        m=metricTemplate(); m.year=year; out=[]; cp=[];
        try
            prefix=fullfile(folder,sprintf('year_%02d',year));
            annualFile=[prefix '_annual_results.mat'];
            check(year,'year_archive_present',isfile(annualFile),annualFile);
            demand(isfile(annualFile),'Missing annual result archive.');
            saved=load(annualFile,'out'); need(saved,{'out'},'annual archive'); out=saved.out;
            need(out,{'version','initialX','finalX','timeDay','history','parameters', ...
                'mechanics','leakage','numerics','pAll','soilNodes','mesh','annualContinuation'},'year output');
            demand(year<=numel(data.coreCheckpoints),'Year checkpoint registration is missing.');
            [cpFile,relocated]=checkpointPath(data.coreCheckpoints{year},folder,year);
            if relocated
                report.warnings{end+1}=sprintf('Year %d checkpoint resolved in relocated case folder: %s.',year,cpFile);
            end
            check(year,'year_checkpoint_present',isfile(cpFile),cpFile);
            demand(isfile(cpFile),'Missing registered year-end checkpoint.');
            cpData=load(cpFile); need(cpData,{'checkpoint','hist','mechanicalHistory', ...
                'jointMechanicalHistory','strainHistory','jointHydraulicHist','jointFreezeHist', ...
                'leakHistSegment','leakHistJoint','leakHistBody','snaps'},'year checkpoint archive');
            cp=cpData.checkpoint;
            need(cp,{'version','X','time','mechanicalReference','signature','pAll','tSoil','qsets','ctri','interface','parameters'},'checkpoint');
            need(cp.mechanicalReference,{'U','T'},'checkpoint installed reference');
            p=out.parameters; mech=out.mechanics; leak=out.leakage;
            if isfield(out,'jointCycleDamage')&&out.jointCycleDamage.enabled
                cyc=out.jointCycleDamage;
                need(cp,{'jointCycleMemory','jointCycleAppliedSeverity','runInventory'},'cycle-damage checkpoint');
                need(cp.runInventory,{'jointCycleHistory','jointCycleIceHistory','startCycleMemory'}, ...
                    'cycle-damage accepted-history inventory');
                check(year,'cycle_checkpoint_full_internal_state', ...
                    isequaln(cp.jointCycleMemory,cyc.finalMemory) ...
                    &&isequal(cp.jointCycleAppliedSeverity,cyc.appliedSeverity), ...
                    'Checkpoint retains the exact applied geometry and complete pending/partial-cycle memory.');
                check(year,'cycle_checkpoint_accepted_history_exact', ...
                    isequal(cp.runInventory.jointCycleHistory,cyc.history) ...
                    &&isequal(cp.runInventory.jointCycleIceHistory,cyc.sampleIceHistory) ...
                    &&isequaln(cp.runInventory.startCycleMemory,cyc.initialMemory), ...
                    'Checkpoint and annual archive retain every accepted local ice state and original year-start memory.');
                if year>1
                    check(year,'cycle_full_memory_inherited',~isempty(prior) ...
                        &&isfield(prior,'jointCycleMemory') ...
                        &&isequaln(cyc.initialMemory,prior.jointCycleMemory), ...
                        'Completed and partial local cycles persist exactly from the preceding year.');
                end
            end
            need(p,{'day','yearDays','durationDays','thm','solver','liningK','leakage'},'core parameters');
            need(mech,{'referenceU','referenceTemperature','history','jointHistory','strainHistory', ...
                'soilStrainDiagnostics','absoluteSoilStrainDiagnostics','installedSoilStrainDiagnostics', ...
                'maximumFreeResidual_N_per_m','omittedNewChannelVolumeFraction', ...
                'exactOpenSlitGeometricVolume_m2','openingGeometry','currentSoilPorosity','mesh','contact'},'THM mechanics');
            demand(finite(out.initialX)&&finite(out.finalX)&&finite(cp.X)&&finite(cp.time), ...
                'Coupled initial/final/checkpoint states must be finite numeric data.');
            h=asMatrix(out.history,{}, {},'accepted-step history');
            demand(size(h,1)>0&&size(h,2)>=14,'Accepted-step history must have fourteen columns and at least one step.');
            mh=asMatrix(mech.history,{}, {},'mechanical history');
            sh=asMatrix(mech.strainHistory,{}, {},'strain history');
            jh=asMatrix(mech.jointHistory,{}, {},'joint mechanical history');
            demand(size(mh,2)==6&&size(sh,2)==7&&size(jh,2)==73 ...
                &&size(mh,1)==size(h,1)&&size(sh,1)==size(h,1)&&size(jh,1)==size(h,1), ...
                'Every accepted step requires mechanical, seven-column strain and eight-joint histories.');
            out.history=h;mech.history=mh;mech.strainHistory=sh;mech.jointHistory=jh;out.mechanics=mech;
            m.acceptedSteps=size(h,1); m.unknowns=numel(out.finalX); m.duration_day=sum(h(:,2));
            if isfield(out,'solverWork')&&isfield(out.solverWork,'rejectedSteps')
                m.rejectedSteps=out.solverWork.rejectedSteps;
                check(year,'solver_work_count',out.solverWork.acceptedSteps==m.acceptedSteps ...
                    &&finite(m.rejectedSteps)&&m.rejectedSteps>=0&&m.rejectedSteps==fix(m.rejectedSteps), ...
                    'Accepted/rejected counts must agree with the stored accepted history.');
            else
                check(year,'solver_work_count',false,'Missing accepted/rejected solver work inventory.');
            end
            check(year,'one_engine_label',strcmp(out.version,report.engineVersion)&&strcmp(cp.version,report.engineVersion), ...
                'Output and checkpoint labels must match the registered annual engine.');
            check(year,'full_365_day_history',p.day==86400&&p.yearDays==365&&p.durationDays==365 ...
                &&all(h(:,2)>0)&&all(diff(h(:,1))>0) ...
                &&near(sum(h(:,2)),365,o.timeTolerance_day,0) ...
                &&near(h(:,1)-[365*(year-1);h(1:end-1,1)],h(:,2),o.timeTolerance_day,0) ...
                &&near(out.timeDay,365*year,o.timeTolerance_day,0) ...
                &&near(cp.time/86400,365*year,o.timeTolerance_day,0), ...
                sprintf('%.12g days, final day %.12g; every interval must be unique and contiguous.',m.duration_day,out.timeDay));
            check(year,'history_time_alignment',isequal(h(:,1),mh(:,1))&&isequal(h(:,1),sh(:,1))&&isequal(h(:,1),jh(:,1)), ...
                'Water/heat, mechanics, joint and strain records use the same accepted steps.');
            check(year,'checkpoint_history_exact',isequal(asMatrix(cpData.hist,{}, {},'checkpoint TH history'),h) ...
                &&isequal(asMatrix(cpData.mechanicalHistory,{}, {},'checkpoint mechanics history'),mh) ...
                &&isequal(asMatrix(cpData.strainHistory,{}, {},'checkpoint strain history'),sh) ...
                &&isequal(asMatrix(cpData.jointMechanicalHistory,{}, {},'checkpoint joint history'),jh), ...
                'Year-end checkpoint preserves all accepted histories exactly.');
            check(year,'final_coupled_state_exact',isequal(out.finalX,cp.X), ...
                'Checkpoint X must equal the year-end T/h/soil-u/ring-u state exactly.');
            if data.schemaVersion>=11
                need(cp,{'runInventory'},'schema-11 complete-year checkpoint');
                need(cp.runInventory,{'startX','startTime','endTime','maxLocalWaterDefect','maxLocalEnergyDefect'},'accepted run inventory');
                check(year,'original_run_inventory_exact',isequal(cp.runInventory.startX,out.initialX) ...
                    &&near(cp.runInventory.startTime/86400,365*(year-1),o.timeTolerance_day,0) ...
                    &&near(cp.runInventory.endTime/86400,365*year,o.timeTolerance_day,0), ...
                    'Original year start X/time and final target remain intact after any month restart.');
            end
            check(year,'installed_reference_exact',isequal(cp.mechanicalReference,data.installedMechanicalReference) ...
                &&isequal(cp.mechanicalReference.U,mech.referenceU) ...
                &&isequal(cp.mechanicalReference.T,mech.referenceTemperature), ...
                'Installed U/T reference is compared exactly, without numerical tolerance or annual reset.');
            need(out.annualContinuation,{'schemaVersion','year','installedReferencePreserved','partialYearResumeUsed'},'annual continuation metadata');
            check(year,'annual_continuation_metadata',out.annualContinuation.schemaVersion==data.schemaVersion ...
                &&out.annualContinuation.year==year&&out.annualContinuation.installedReferencePreserved ...
                &&(~out.annualContinuation.partialYearResumeUsed||data.schemaVersion>=11), ...
                'A schema-11 month restart qualifies only when its merged accepted history and original start X pass all full-year checks.');
            check(year,'checkpoint_parameter_exact',isequaln(cp.parameters,p), ...
                'Output and checkpoint core parameters must agree.');
            geometry=struct('pAll',cp.pAll,'tSoil',cp.tSoil,'qsets',{cp.qsets}, ...
                'ctri',cp.ctri,'interface',cp.interface,'soilNodes',out.soilNodes);
            if year==1
                firstSignature=cp.signature; firstGeometry=geometry;
                report.parameters=p; report.parameterSummary=parameterSummary(p,out);
                N=size(out.pAll,1); ns=numel(out.soilNodes);
                check(year,'installation_state_exact',isequal(out.initialX(1:N),cp.mechanicalReference.T) ...
                    &&isequal(out.initialX(N+ns+1:end),cp.mechanicalReference.U), ...
                    'Fresh initial coupled state contains the actual installed T/U equilibrium.');
            else
                check(year,'prior_entire_state_exact',~isempty(prior)&&isequal(out.initialX,prior.X), ...
                    'Every new year inherits the entire previous year-end X exactly.');
                check(year,'physical_parameters_unchanged',isequaln(cp.signature,firstSignature), ...
                    'Checkpoint physical signatures remain identical across all years.');
                check(year,'mesh_unchanged',isequaln(geometry,firstGeometry), ...
                    'Thermal/hydraulic and soil mechanical geometry remain identical across years.');
            end
            check(year,'enabled_two_way_THM',p.thm.enabled&&~strcmp(p.initial.mode,'spinup'), ...
                'A TH-only periodic shortcut does not qualify as a full THM history.');
            check(year,'eight_independent_mechanical_blocks',eightBlocks(mech.mesh), ...
                'Eight disconnected Q4 segment blocks, independent coincident cut nodes, and no element crossing a cut.');

            [leakChecks,volumes,potential]=inspectLeakage(out,year,o);
            append(year,leakChecks);
            [diffusionChecks,diffusion]=inspectDiffusion(out,year);
            append(year,diffusionChecks);
            hydraulicAudits=[hydraulicAudits;diffusion]; %#ok<AGROW>
            check(year,'checkpoint_hydraulics_exact',isequal(cpData.jointHydraulicHist,leak.hydraulics.history), ...
                'Checkpoint preserves the complete accepted hydraulic history exactly, including resumed prefixes.');
            check(year,'checkpoint_joint_leak_freeze_exact',isequal(cpData.leakHistSegment,leak.segmentRate_m2_per_s) ...
                &&isequal(cpData.leakHistJoint,leak.jointRate_m2_per_s) ...
                &&isequal(cpData.leakHistBody,leak.bodyRate_m2_per_s) ...
                &&isfield(leak,'jointFreeze')&&isequal(cpData.jointFreezeHist,leak.jointFreeze.history), ...
                'Checkpoint preserves every leakage/freezing accepted prefix exactly.');
            m.totalLeak_m3_per_m=volumes.total; m.jointLeak_m3_per_m=sum(volumes.joint);
            m.bodyLeak_m3_per_m=volumes.body; m.sameStateFarHeadLeak_m3_per_m=sum(potential.annual);
            m.supplyDifference_m3_per_m=m.sameStateFarHeadLeak_m3_per_m-m.totalLeak_m3_per_m;
            m.supplyLossFraction=divideDefined(m.supplyDifference_m3_per_m,m.sameStateFarHeadLeak_m3_per_m);
            m.maxOpenSlitConductance_m_s=max(jh(:,58:65),[],'all');
            m.maxCapillaryConductance_m_s=max(jh(:,50:57),[],'all');
            m.openSlitActiveDuration_day=sum(h(any(jh(:,58:65)>0,2),2));
            share=zeros(size(jh,1),8);active=jh(:,66:73)>0;
            open=jh(:,58:65);totalC=jh(:,66:73);share(active)=open(active)./totalC(active);
            m.maxOpenConductanceShare=max(share,[],'all');
            aggregateChecks=inspectAggregate(ledger,leak,year,o);
            append(year,aggregateChecks);

            need(out.numerics,{'maximumAcceptedLocalWaterDefect_m2','maximumAcceptedLocalEnergyDefect_Jm'},'local water/heat ledgers');
            [soilV,thermalV,q4Orientation]=nodalVolumes(out);
            m.positiveQ4Elements=q4Orientation.positiveElements;
            m.negativeQ4Elements=q4Orientation.negativeElements;
            check(year,'consistent_nonsingular_Q4_orientation',true, ...
                sprintf('Each element has four nonzero, same-sign Gauss Jacobians; %d positive and %d negative orientations use abs(detJ).',m.positiveQ4Elements,m.negativeQ4Elements));
            m.waterSignedAnnualDefect_m2=sum(h(:,8)); m.energySignedAnnualDefect_J_per_m=sum(h(:,9));
            m.waterAbsoluteAnnualDefect_m2=sum(abs(h(:,8))); m.energyAbsoluteAnnualDefect_J_per_m=sum(abs(h(:,9)));
            m.maxGlobalWaterDefect_m2=max(abs(h(:,8))); m.maxGlobalEnergyDefect_J_per_m=max(abs(h(:,9)));
            m.maxLocalWaterDefect_m2=out.numerics.maximumAcceptedLocalWaterDefect_m2;
            m.maxLocalEnergyDefect_J_per_m=out.numerics.maximumAcceptedLocalEnergyDefect_Jm;
            m.maxScaledResidual=max(h(:,11)); m.maxMechanicalResidual_N_per_m=max(mh(:,5));
            % Accepted scaled infinity residual bounds each nodal ledger by
            % tol*its nodal volume*scale. Global bounds sum these volumes.
            m.globalWaterTolerance_m2=p.solver.tol*p.solver.massScale*sum(soilV);
            m.localWaterTolerance_m2=p.solver.tol*p.solver.massScale*max(soilV);
            m.globalEnergyTolerance_J_per_m=p.solver.tol*p.solver.energyScale*sum(thermalV);
            m.localEnergyTolerance_J_per_m=p.solver.tol*p.solver.energyScale*max(thermalV);
            m.mechanicalTolerance_N_per_m=p.solver.tol/p.thm.forceResidualScale;
            waterRound=128*eps(max(1,max(abs(h(:,4)))));
            energyRound=128*eps(max(1,max(abs(h(:,5)))));
            check(year,'scaled_newton_residual',all(h(:,11)>=0)&&m.maxScaledResidual<=p.solver.tol*(1+1e-7), ...
                sprintf('max %.6g, solver tol %.6g.',m.maxScaledResidual,p.solver.tol));
            check(year,'global_water_ledger',m.maxGlobalWaterDefect_m2<=m.globalWaterTolerance_m2+waterRound, ...
                sprintf('max %.6g m2; solver-volume bound %.6g m2 (+ roundoff %.3g).',m.maxGlobalWaterDefect_m2,m.globalWaterTolerance_m2,waterRound));
            check(year,'global_heat_ledger',m.maxGlobalEnergyDefect_J_per_m<=m.globalEnergyTolerance_J_per_m+energyRound, ...
                sprintf('max %.6g J/m; solver-volume bound %.6g J/m (+ roundoff %.3g).',m.maxGlobalEnergyDefect_J_per_m,m.globalEnergyTolerance_J_per_m,energyRound));
            check(year,'local_water_ledger',finite(m.maxLocalWaterDefect_m2)&&m.maxLocalWaterDefect_m2>=0 ...
                &&m.maxLocalWaterDefect_m2<=m.localWaterTolerance_m2+waterRound, ...
                sprintf('max %.6g m2; largest-node bound %.6g m2.',m.maxLocalWaterDefect_m2,m.localWaterTolerance_m2));
            check(year,'local_heat_ledger',finite(m.maxLocalEnergyDefect_J_per_m)&&m.maxLocalEnergyDefect_J_per_m>=0 ...
                &&m.maxLocalEnergyDefect_J_per_m<=m.localEnergyTolerance_J_per_m+energyRound, ...
                sprintf('max %.6g J/m; largest-node bound %.6g J/m.',m.maxLocalEnergyDefect_J_per_m,m.localEnergyTolerance_J_per_m));
            check(year,'mechanical_equilibrium',all(mh(:,5)>=0) ...
                &&m.maxMechanicalResidual_N_per_m<=m.mechanicalTolerance_N_per_m*(1+1e-6) ...
                &&near(mech.maximumFreeResidual_N_per_m,mh(end,5),1e-10,o.relativeTolerance), ...
                sprintf('max free residual %.6g N/m; unscaled solver bound %.6g N/m.',m.maxMechanicalResidual_N_per_m,m.mechanicalTolerance_N_per_m));
            m.maxIncrementalVolumetricStrain=max(sh(:,2)); m.maxIncrementalPrincipalStrain=max(sh(:,3));
            m.maxIncrementalTensorShearStrain=max(sh(:,4)); m.maxIncrementalGradientNorm=max(sh(:,5));
            m.maxAbsolutePrincipalStrain=max(sh(:,6)); m.maxAbsoluteGradientNorm=max(sh(:,7));
            m.installedPrincipalStrain=mech.installedSoilStrainDiagnostics.maxAbsPrincipalStrain;
            m.installedGradientNorm=mech.installedSoilStrainDiagnostics.maxGradientNorm;
            m.maxOmittedSlitVolumeFraction=max(mh(:,6));
            m.minFinalPorosity=min(mech.currentSoilPorosity); m.maxFinalPorosity=max(mech.currentSoilPorosity);
            check(year,'full_strain_history_domain',all(sh(:,2:end)>=0,'all') ...
                &&m.maxIncrementalVolumetricStrain<=p.thm.maximumVolumetricStrain ...
                &&m.maxIncrementalPrincipalStrain<=p.thm.maximumPrincipalStrain ...
                &&m.maxAbsolutePrincipalStrain<=p.thm.maximumPrincipalStrain ...
                &&m.maxIncrementalGradientNorm<=p.thm.maximumDisplacementGradient ...
                &&m.maxAbsoluteGradientNorm<=p.thm.maximumDisplacementGradient, ...
                'Every accepted step is checked for incremental volume/principal/shear/gradient and absolute principal/gradient values.');
            check(year,'installed_absolute_strain_domain',finite(m.installedPrincipalStrain)&&finite(m.installedGradientNorm) ...
                &&m.installedPrincipalStrain>=0&&m.installedGradientNorm>=0 ...
                &&m.installedPrincipalStrain<=p.thm.maximumPrincipalStrain ...
                &&m.installedGradientNorm<=p.thm.maximumDisplacementGradient, ...
                'Absolute installation deformation also satisfies the reference-geometry validity guards.');
            check(year,'final_strain_consistency',near(mh(:,2),sh(:,2),1e-10,o.relativeTolerance) ...
                &&near(sh(end,2:5),[mech.soilStrainDiagnostics.maxAbsVolumetricStrain, ...
                mech.soilStrainDiagnostics.maxAbsPrincipalStrain,mech.soilStrainDiagnostics.maximumTensorShearStrain, ...
                mech.soilStrainDiagnostics.maxGradientNorm],1e-10,o.relativeTolerance) ...
                &&near(sh(end,6:7),[mech.absoluteSoilStrainDiagnostics.maxAbsPrincipalStrain, ...
                mech.absoluteSoilStrainDiagnostics.maxGradientNorm],1e-10,o.relativeTolerance), ...
                'Final full-strain diagnostics agree with the accepted history.');
            check(year,'positive_physical_porosity',finite(mech.currentSoilPorosity) ...
                &&m.minFinalPorosity>0&&m.maxFinalPorosity<1, ...
                sprintf('Final nodal porosity range %.6g..%.6g.',m.minFinalPorosity,m.maxFinalPorosity));
            check(year,'slit_inventory_history_domain',all(mh(:,6)>=0) ...
                &&m.maxOmittedSlitVolumeFraction<=p.thm.maxOmittedSlitVolumeFraction ...
                &&near(mech.omittedNewChannelVolumeFraction,mh(end,6),1e-12,o.relativeTolerance), ...
                sprintf('Largest single-instant exact geometric inventory fraction %.6g; bound %.6g. This is not a leakage-error bound.',m.maxOmittedSlitVolumeFraction,p.thm.maxOmittedSlitVolumeFraction));
            [exactV,baseV]=independentSlitVolume(mech,p);
            frac=zeros(8,1); pos=baseV>0; frac(pos)=exactV(pos)./baseV(pos);
            check(year,'exact_positive_part_slit_geometry',near(exactV,mech.exactOpenSlitGeometricVolume_m2,1e-16,o.relativeTolerance) ...
                &&all(exactV(~pos)==0)&&near(max(frac),mh(end,6),1e-12,o.relativeTolerance), ...
                'Independent nodal P1 positive-part radial integration; no temperature-midpoint/minimum-gap proxy.');
            if isfield(out,'snapshots')&&~isempty(out.snapshots)
                validSnapshots=true;
                for s=1:numel(out.snapshots)
                    validSnapshots=validSnapshots&&isfield(out.snapshots(s),'U') ...
                        &&finite(out.snapshots(s).U)&&numel(out.snapshots(s).U)==numel(mech.referenceU);
                end
                check(year,'snapshot_full_mechanical_state',validSnapshots, ...
                    'Snapshots retain soil and ring mechanical DOFs, rather than only diagnostic scalar strain.');
            else
                check(year,'snapshot_full_mechanical_state',false,'Missing full mechanical snapshots.');
            end
            if year==20
                report.year20=year20Summary(out,potential,exactV,baseV);
                lastOut=out;
            end
        catch ME
            check(year,'readable_complete_year_data',false,sprintf('%s: %s',ME.identifier,ME.message));
        end
        metrics(end+1,1)=m; %#ok<AGROW>
        % Never connect a later year to an unreadable or unregistered state.
        if isstruct(cp)&&isfield(cp,'X')
            prior=struct('X',cp.X);
            if isfield(cp,'jointCycleMemory'), prior.jointCycleMemory=cp.jointCycleMemory; end
        else
            prior=[];
        end
    end
    if n>0
        check(0,'registered_last_year_exact',isstruct(data.lastYear)&&isfield(data.lastYear,'finalX') ...
            &&~isempty(prior)&&isequal(data.lastYear.finalX,prior.X), ...
            'Aggregate lastYear must agree with the final registered year-end checkpoint.');
    end
    if n==20 && ~isempty(lastOut)
        check(0,'last_year_7300_day',near(lastOut.timeDay,7300,o.timeTolerance_day,0), ...
            'The twentieth-year response follows 7300 actual simulated days.');
        for key={'year20AnnualLeakage_L_per_m','year20ByJoint_L_per_m','year20ByMonth_L_per_m'}
            if isfield(data,key{1})
                target=report.year20.(key{1});
                check(0,['published_' key{1}],finite(data.(key{1})) ...
                    &&near(data.(key{1})(:),target(:),1000*o.volumeAbsTolerance_m2,o.relativeTolerance), ...
                    'Published year-20 summaries match the original conservative accepted-step history.');
            end
        end
    end
catch ME
    check(0,'readable_case_data',false,sprintf('%s: %s',ME.identifier,ME.message));
end
report.yearSummary=struct2table(metrics);
report.hydraulicAudit=hydraulicAudits;
report.checks=struct2table(checks);
if ~isempty(checks)
    failed=find(~[checks.passed]);
    for k=failed
        report.issues{end+1}=sprintf('Year %d / %s: %s',checks(k).year,checks(k).name,checks(k).detail);
    end
end
report.passed=isequal(report.completedYears,20)&&numel(metrics)==20 ...
    &&~isempty(checks)&&all([checks.passed])&&~isempty(lastOut);
if report.passed&&o.makePlots
    try
        figs=plotReport(report,o);
    catch ME
        report.warnings{end+1}=sprintf('Plot export unavailable: %s',ME.message);
    end
elseif o.makePlots
    report.plotSuppressionReason='Scientific year-20 plots require a complete passing 20-year numerical audit.';
end
if o.saveReport
    destination=char(o.outputFolder);
    if isempty(destination)&&~isempty(folder), destination=fullfile(folder,'validation_20year'); end
    if isempty(destination)
        report.warnings{end+1}='Report save skipped because no output folder is available.';
    else
        try
            if ~isfolder(destination), mkdir(destination); end
            exportReport(report,figs,destination);
        catch ME
            report.warnings{end+1}=sprintf('Report save unavailable: %s',ME.message);
        end
    end
end
fprintf('THM 20-year audit: passed=%d; completed=%g/20; failed checks=%d.\n', ...
    report.passed,report.completedYears,numel(report.issues));

    function check(year,name,passed,detail)
        checks(end+1)=struct('year',year,'name',name,'passed',logical(passed),'detail',char(detail));
    end
    function append(year,list)
        for i=1:numel(list),check(year,list(i).name,list(i).passed,list(i).detail);end
    end
end

function [data,folder,file]=readInput(input)
file=''; folder='';
if isstruct(input)&&isscalar(input)
    data=input;
    if isfield(data,'config')&&isfield(data.config,'outputFolder'), folder=char(data.config.outputFolder); end
elseif (ischar(input)&&isrow(input))||(isstring(input)&&isscalar(input))
    path=char(input);
    path=canonicalPath(path);
    if isfolder(path)
        folder=path;
        % Progress is the live authoritative year-commit inventory. A final
        % result may be older after a resumed run is interrupted mid-year.
        file=fullfile(folder,'thm_progress.mat');
        if ~isfile(file),file=fullfile(folder,'TH_segment_year20_thm_results.mat');end
    else
        file=path; folder=fileparts(file);
    end
    demand(isfile(file),sprintf('No annual progress/result MAT found: %s.',file));
    saved=load(file);
    if isfield(saved,'progress'),data=saved.progress;
    elseif isfield(saved,'result'),data=saved.result;
    else,error('THM20Audit:Input','MAT must contain annual progress or result.');end
else
    error('THM20Audit:Input','Supply a scalar annual result/progress struct, output folder, or MAT path.');
end
demand(isstruct(data)&&isscalar(data),'Annual data must be a scalar struct.');
if ~isempty(folder),folder=canonicalPath(folder);end
end

function [file,relocated]=checkpointPath(registered,folder,year)
demand((ischar(registered)&&isrow(registered))||(isstring(registered)&&isscalar(registered)), ...
    'Checkpoint registration must be a path string.');
registered=char(registered); [~,name,extension]=fileparts(registered);
demand(strcmp([name extension],sprintf('year_%02d_checkpoint.mat',year)), ...
    'Registered checkpoint filename does not identify this complete year.');
file=registered; relocated=false;
local=fullfile(folder,[name extension]);
% An explicitly supplied case folder is authoritative after relocation.
if isfile(local)
    file=local;
    relocated=~strcmp(canonicalPath(local),canonicalPath(registered));
end
end

function path=canonicalPath(path)
% Java's working directory need not track MATLAB cd/pwd. Anchor relative
% inputs to the calling MATLAB directory before invoking Java canonicalize.
path=char(path);
if ~java.io.File(path).isAbsolute(),path=fullfile(pwd,path);end
path=char(java.io.File(path).getCanonicalPath());
end

function expected=columnSchema()
expected.annual=struct('metadata','annualColumns','columns',{{'year','totalLeak_m3_per_m', ...
    'jointLeak_m3_per_m','bodyLeak_m3_per_m','waterDefect_m3_per_m','energyDefect_J_per_m', ...
    'maxAbsWaterDefect_m3_per_m','maxAbsEnergyDefect_J_per_m','acceptedSteps','maxScaledResidual', ...
    'maxAbsVolumetricStrain','maxJointGap_m','maxFreeMechanicalResidual_N_per_m','maxOmittedNewSlitVolumeFraction'}});
jointNames=arrayfun(@(j)sprintf('joint%dLeak_m3_per_m',j),1:8,'UniformOutput',false);
expected.monthly=struct('metadata','monthlyColumns','columns',{{'year','month','totalLeak_m3_per_m',jointNames{:}}});
expected.byJoint=struct('metadata','jointColumns','columns',{{'year','jointID','angle_rad','leak_m3_per_m'}});
expected.daily=struct('metadata','dailyColumns','columns',{{'year','dayOfYear','absoluteEndDay','totalLeak_m3_per_m',jointNames{:}}});
end

function A=asMatrix(value,columns,required,label)
if istable(value)
    columns=value.Properties.VariableNames; A=table2array(value);
else
    A=value;
end
demand(isnumeric(A)&&isreal(A)&&ismatrix(A)&&all(isfinite(A(:))), ...
    sprintf('%s contains missing/nonfinite/nonnumeric values.',label));
A=double(A);
if ~isempty(required)
    columns=cellstr(string(columns));
    demand(numel(columns)==size(A,2),sprintf('%s column metadata does not match the matrix.',label));
    [present,index]=ismember(required,columns);
    demand(all(present)&&numel(unique(columns))==numel(columns),sprintf('%s is missing required unique column names.',label));
    A=A(:,index);
end
end

function [checks,volume,potential]=inspectLeakage(out,year,o)
checks=emptyChecks(); h=out.history; leak=out.leakage; p=out.parameters; n=size(h,1);
need(leak,{'segmentRate_m2_per_s','jointRate_m2_per_s','bodyRate_m2_per_s','totalVolume_m3_per_m', ...
    'volumeByJoint_m3_per_m','volumeBySegment_m3_per_m','bodyVolume_m3_per_m', ...
    'dailyBinEndDay','dailyTotalVolume_m3_per_m','dailyVolumeByJoint_m3_per_m','dailyBodyVolume_m3_per_m', ...
    'monthlyTotalVolume_m3_per_m','monthlyVolumeByJoint_m3_per_m','jointAnglesRad','hydraulics','jointQuadrature'},'leakage ledger');
seg=asMatrix(leak.segmentRate_m2_per_s,{}, {},'segment rates');
joint=asMatrix(leak.jointRate_m2_per_s,{}, {},'joint rates');
body=asMatrix(leak.bodyRate_m2_per_s,{}, {},'body rates');
demand(isequal(size(seg),[n 8])&&isequal(size(joint),[n 8])&&isequal(size(body),[n 1]), ...
    'Rates must cover every accepted step and all eight joints/sectors.');
dt=h(:,2)*86400; t=[365*(year-1);h(:,1)];
volume=struct('total',sum(h(:,10).*dt),'joint',sum(joint.*dt,1),'body',sum(body.*dt));
checks(end+1)=item('nonnegative_outflow',all(seg(:)>=-1e-14)&&all(joint(:)>=-1e-14)&&all(body(:)>=-1e-14), ...
    'All accepted sector/joint/body outflows must be nonnegative.');
checks(end+1)=item('step_leak_partition',near(sum(seg,2),h(:,10),o.rateAbsTolerance_m2_s,o.relativeTolerance) ...
    &&near(sum(joint,2)+body,h(:,10),o.rateAbsTolerance_m2_s,o.relativeTolerance), ...
    'Every accepted rate closes both sector and joint/body partitions.');
checks(end+1)=item('integrated_annual_leak',near(volume.total,leak.totalVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(volume.joint(:),leak.volumeByJoint_m3_per_m(:),o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(sum(seg.*dt,1)',leak.volumeBySegment_m3_per_m(:),o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(volume.body,leak.bodyVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance), ...
    'Annual totals are recomputed from original accepted dt, in m3 per metre of tunnel axis.');
daily=acceptedBins(t,joint,365*(year-1)+(0:365)');
dailyBody=acceptedBins(t,body,365*(year-1)+(0:365)');
monthEdges=365*(year-1)+[0 cumsum([31 28 31 30 31 30 31 31 30 31 30 31])];
monthly=acceptedBins(t,joint,monthEdges');
monthlyTotal=sum(monthly,2)+acceptedBins(t,body,monthEdges');
checks(end+1)=item('daily_joint_closure',near(daily,leak.dailyVolumeByJoint_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(dailyBody,leak.dailyBodyVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(sum(daily,2)+dailyBody,leak.dailyTotalVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(leak.dailyBinEndDay,365*(year-1)+(1:365)',o.timeTolerance_day,0), ...
    'All 365 daily bins close using conservative overlap integration of accepted intervals.');
checks(end+1)=item('monthly_joint_closure',near(monthly,leak.monthlyVolumeByJoint_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(monthlyTotal,leak.monthlyTotalVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(sum(monthlyTotal),volume.total,o.volumeAbsTolerance_m2,o.relativeTolerance), ...
    'All twelve month and eight-joint bins close with the repeated 365-day calendar.');
jointLaw=joint_microchannel_damage(p.leakage);
Cref=expand(jointLaw.jointConductance,8); sealed=Cref==0;
CrefHistory=repmat(Cref',n,1);
if isfield(out,'jointCycleDamage')&&out.jointCycleDamage.enabled
    cycleAudit=audit_joint_cycle_damage(out);
    checks(end+1)=item('accepted_local_cycle_memory_and_geometry',cycleAudit.passed, ...
        'Saved sample ice histories reconstruct complete local cycles, queued/applied severity and constant connected volume.');
    CrefHistory=out.jointCycleDamage.history(:,42:49);
end
checks(end+1)=item('impermeable_concrete_body',p.liningK==0&&all(body==0)&&volume.body==0&&leak.bodyVolume_m3_per_m==0, ...
    'Concrete body must have K=0 and exactly zero leakage.');
checks(end+1)=item('sealed_joints_do_not_create_leaks',all(joint(:,sealed)==0,'all'), ...
    'Cref=0 joints remain sealed; mechanical opening cannot invent an uncalibrated gasket failure.');
hyd=leak.hydraulics;
need(hyd,{'history','drivingHead_m','referenceHead_m','unrestrictedRate_m2_per_s', ...
    'waterTableDepth_m','minTrans_m_per_s','negativeTransEdges'},'hydraulic accepted history');
raw=asMatrix(hyd.history,{}, {},'hydraulic history');
qr=asMatrix(hyd.unrestrictedRate_m2_per_s,{}, {},'same-state far-head rate');
demand(isequal(size(raw),[n 27])&&isequal(size(qr),[n 8]),'Hydraulic history must contain three eight-joint blocks plus three diagnostics.');
checks(end+1)=item('hydraulic_history_exact',isequal(raw(:,1:8),hyd.drivingHead_m) ...
    &&isequal(raw(:,9:16),hyd.referenceHead_m)&&isequal(raw(:,17:24),qr) ...
    &&isequal(raw(:,25),hyd.waterTableDepth_m)&&isequal(raw(:,26),hyd.minTrans_m_per_s) ...
    &&isequal(raw(:,27),hyd.negativeTransEdges), 'Named hydraulic fields preserve the original record exactly.');
checks(end+1)=item('signed_FEM_edge_diagnostics',all(raw(:,27)>=0&raw(:,27)==fix(raw(:,27)))&&all(qr(:)>=0), ...
    'Signed aggregate P1 edge weights are retained, not clipped. Nonnegative element K and diffusive-energy/physical-subface identities are checked separately; M-matrix/DMP is not guaranteed.');
Ctotal=out.mechanics.jointHistory(:,66:73); Cbase=out.mechanics.jointHistory(:,50:57);
Copen=out.mechanics.jointHistory(:,58:65);
need(leak,{'jointFreeze'},'accepted freeze/viscosity history');
need(leak.jointFreeze,{'history','factor','minTemperature_C','maxTemperature_C','frozenLengthFraction'},'accepted freeze/viscosity history');
freeze=asMatrix(leak.jointFreeze.history,{}, {},'freeze/viscosity history');
demand(isequal(size(freeze),[n 32]),'Freezing history must contain four eight-joint accepted blocks.');
checks(end+1)=item('same_state_capillary_phase_conductance',all(freeze(:,1:8)>=0,'all') ...
    &&all(freeze(:,9:16)<=freeze(:,17:24),'all') ...
    &&all(freeze(:,25:32)>=0&freeze(:,25:32)<=1,'all') ...
    &&isequal(freeze(:,1:8),leak.jointFreeze.factor) ...
    &&near(freeze(:,1:8).*CrefHistory,Cbase,1e-20,o.relativeTolerance), ...
    'Accepted capillary conductance uses the same recorded thermal freeze/viscosity factor; open flow is a separate physical branch.');
checks(end+1)=item('total_same_state_conductance',all(Cbase(:)>=0)&&all(Copen(:)>=0)&&all(Ctotal(:)>=0) ...
    &&near(Cbase+Copen,Ctotal,1e-20,o.relativeTolerance)&&all(Ctotal(:,sealed)==0,'all'), ...
    'Ceff is the recorded capillary-plus-open-slit conductance at the same accepted state.');
Q=leak.jointQuadrature; need(Q,{'phi','weight','jointID'},'Robin quadrature');
phi=Q.phi; w=double(Q.weight(:)); id=double(Q.jointID(:)); nq=numel(w);
demand(size(phi,1)==nq&&size(phi,2)==numel(out.soilNodes)&&all(isfinite(nonzeros(phi))) ...
    &&all(w>0)&&all(isfinite(w))&&numel(id)==nq&&all(id>=1&id<=8&id==fix(id)), ...
    'Invalid stored joint Robin quadrature.');
J=sparse((1:nq)',id,1,nq,8); xy=out.pAll(out.soilNodes,:); center=out.mechanics.mesh.center;
radius=hypot(xy(:,1)-center(1),xy(:,2)-center(2)); drop=zeros(size(radius)); nonzero=radius>0;
drop(nonzero)=-p.liningThickness*(xy(nonzero,2)-center(2))./radius(nonzero);
constant=phi*(xy(:,2)-min(out.pAll(:,2))+drop-p.leakage.interiorHead);
head=constant'-hyd.waterTableDepth_m; entry=expand(p.leakage.jointEntryHead,8);
excess=max(head-entry(id)',0);
reconstructed=(Ctotal(:,id).*excess.*w')*J;
checks(end+1)=item('same_state_far_head_reference',near(head.*w'*J,hyd.referenceHead_m,1e-9,o.relativeTolerance) ...
    &&near(reconstructed,qr,1e-20,o.relativeTolerance), ...
    'Recomputed Ceff times FAR head with entry clipping inside each actual Robin quadrature point; local supply drawdown is excluded.');
potential=struct('annual',sum(qr.*dt,1),'daily',acceptedBins(t,qr,365*(year-1)+(0:365)'), ...
    'monthly',acceptedBins(t,qr,monthEdges'),'Ceff_m_s',Ctotal,'rate_m2_s',qr, ...
    'referenceHead_m',hyd.referenceHead_m,'waterTableDepth_m',hyd.waterTableDepth_m);
end

function [checks,d]=inspectDiffusion(out,year)
% P1 diffusion is an element-energy operator, not necessarily an M-matrix.
% Verify its factorization and physical median-dual divergence independently
% of the signed cotangent decomposition; never clip negative edge weights.
need(out,{'state','hydro'},'year-end constitutive/physical flux state');
need(out.state,{'elementK','quadrature'},'year-end constitutive state');
need(out.state.quadrature,{'K'},'SCV permeability state');
need(out.hydro,{'subfaceFlux_m2s'},'physical median-dual face flux');
xy=out.pAll;tri=out.mesh.tSoil;nodes=out.soilNodes(:);ns=numel(nodes);ne=size(tri,1);
Ke=out.state.elementK(:);qK=out.state.quadrature.K;
demand(finite(Ke)&&numel(Ke)==ne&&finite(qK),'Year-end element/SCV K must be finite.');
cv=build_soil_cvfe_flux(xy,tri,nodes);
P=build_soil_scv_quadrature(xy(nodes,:),cv.localTri, ...
    out.parameters.integration.soilQuadratureOrder,out.parameters.integration.subdivisionLevel);
averageK=P.elementAverage*qK(:);
% Coordinate-difference formula provides a second element gradient assembly,
% while the public CVFE kernel independently uses direct nodal coordinates.
a=xy(tri(:,1),:);b=xy(tri(:,2),:)-a;c=xy(tri(:,3),:)-a;
signed=b(:,1).*c(:,2)-b(:,2).*c(:,1);area=abs(signed)/2;
gx=[b(:,2)-c(:,2),c(:,2),-b(:,2)]./signed;
gz=[c(:,1)-b(:,1),-c(:,1),b(:,1)]./signed;
local=cv.localTri;ii=zeros(ne,9);jj=ii;values=ii;slot=0;
for i=1:3
    for j=1:3
        slot=slot+1;ii(:,slot)=local(:,i);jj(:,slot)=local(:,j);
        values(:,slot)=area.*Ke.*(gx(:,i).*gx(:,j)+gz(:,i).*gz(:,j));
    end
end
K=sparse(ii(:),jj(:),values(:),ns,ns);
rowDegree=full(max(sum(spones(K),2)));normK=norm(K,inf);
coefficientRound=512*eps*max(1,rowDegree)*max(normK,realmin);
checks=emptyChecks();
checks(end+1)=item('nonnegative_element_and_SCV_permeability',all(Ke>=0)&&all(qK(:)>=0) ...
    &&norm(Ke-averageK,inf)<=128*eps*max(norm(Ke,inf),realmin), ...
    'Year-end element K is a positive SCV average of nonnegative physical permeability; signed nodal edge coefficients are a separate quantity.');
% sqrt weights are used only after a nonnegative-K check. A negative raw K
% fails rather than being clipped or silently turned into a complex factor.
demand(all(Ke>=0),'Negative element permeability prevents a real diffusive-energy factorization.');
element=repmat((1:ne)',3,1);
Gx=sparse(element,local(:),gx(:),ne,ns);Gz=sparse(element,local(:),gz(:),ne,ns);
W=spdiags(sqrt(area.*Ke),0,ne,ne);L=[W*Gx;W*Gz];
factorError=norm(K-L'*L,inf);symmetryError=norm(K-K',inf);
constantError=norm(K*ones(ns,1),inf);
H=out.h(:)-xy(nodes,2);
gradient=[Gx*H,Gz*H];elementEnergy=sum(area.*Ke.*sum(gradient.^2,2));
matrixEnergy=H'*(K*H);
energyRound=coefficientRound*max(1,norm(H,1)*norm(H,inf));
checks(end+1)=item('diffusive_element_energy_factorization',factorError<=coefficientRound ...
    &&symmetryError<=coefficientRound&&elementEnergy>=0 ...
    &&matrixEnergy>=-energyRound&&abs(matrixEnergy-elementEnergy)<=energyRound, ...
    sprintf('K_FE=L''L, L=sqrt(Ae*Ke)*grad; factor error %.3g, energy %.6g, arithmetic bound %.3g. Coercivity is modulo component constants/Dirichlet anchoring.',factorError,elementEnergy,energyRound));
checks(end+1)=item('constant_hydraulic_potential_patch',constantError<=coefficientRound, ...
    sprintf('K_FE*1 infinity norm %.3g, coefficient-arithmetic bound %.3g; gravity-balanced constant total potential gives zero diffusion.',constantError,coefficientRound));
[F,DH]=soil_cvfe_flux(out.h,Ke,cv,'element');
operatorError=norm(cv.inc*DH-K,inf);divergenceError=norm(cv.inc*F-K*H,inf);
flowRound=coefficientRound*max(1,norm(H,inf));
recorded=out.hydro.subfaceFlux_m2s(:);
checks(end+1)=item('physical_subface_FEM_divergence_identity',finite(recorded)&&numel(recorded)==numel(F) ...
    &&operatorError<=coefficientRound&&divergenceError<=flowRound ...
    &&norm(recorded-F,inf)<=flowRound, ...
    sprintf('Actual Darcy median-dual face flux, inc*D_face=K_FE; operator error %.3g, state divergence error %.3g. No cotangent edge sign is used for heat donors.',operatorError,divergenceError));

pairs=[local(:,[1 2]);local(:,[2 3]);local(:,[3 1])];
base=-[area.*(gx(:,1).*gx(:,2)+gz(:,1).*gz(:,2)); ...
    area.*(gx(:,2).*gx(:,3)+gz(:,2).*gz(:,3)); ...
    area.*(gx(:,3).*gx(:,1)+gz(:,3).*gz(:,1))];
[edge,~,map]=unique(sort(pairs,2),'rows');nedge=size(edge,1);
geometric=accumarray(map,base,[nedge 1]);
trans=accumarray(map,base.*repmat(Ke,3,1),[nedge 1]);
edgeArithmetic=64*eps*accumarray(map,abs(base).*repmat(Ke,3,1),[nedge 1]);
globalEdge=nodes(edge);lining=false(nedge,1);
if isfield(out.mesh,'soilTunnelCoupleTable')&&~isempty(out.mesh.soilTunnelCoupleTable)
    lining=ismember(sort(globalEdge,2),sort(out.mesh.soilTunnelCoupleTable(:,[2 3]),2),'rows');
end
geometricNegative=geometric< -1e-12;finalNegative=trans<0;
hyd=out.leakage.hydraulics;history=hyd.history;
certainNegative=sum(trans<-edgeArithmetic);ambiguous=sum(abs(trans)<=edgeArithmetic);
recordedCount=history(end,27);
checks(end+1)=item('signed_FEM_edge_inventory_consistency',out.numerics.negativeHydraulicEdges==sum(geometricNegative) ...
    &&abs(out.numerics.minimumGeometricEdgeWeight-min(geometric))<=128*eps*max(norm(geometric,inf),realmin) ...
    &&abs(history(end,26)-min(trans))<=max(edgeArithmetic) ...
    &&recordedCount>=certainNegative&&recordedCount<=certainNegative+ambiguous, ...
    sprintf('Geometry negative edges %d; final strictly negative %d (recorded %g); only cancellation signs within per-edge roundoff are ambiguous.',sum(geometricNegative),sum(finalNegative),recordedCount));
checks(end+1)=item('lining_negative_weight_protection',out.numerics.negativeLiningEdges==0 ...
    &&~any(geometricNegative&lining)&&all(trans(lining)>=-edgeArithmetic(lining)), ...
    sprintf('Lining-support negative geometric count %d; final materially negative count %d. Robin support/weights remain positive.',sum(geometricNegative&lining),sum(lining&trans<-edgeArithmetic)));
liningContribution=lining(map);
allLiningElementTermsNonnegative=all(base(liningContribution)>=0);
negativeLiningTerms=liningContribution&base<0;
termElement=repmat((1:ne)',3,1);termGlobalEdge=nodes(sort(pairs,2));
liningNegativeContributions=table(termElement(negativeLiningTerms),termGlobalEdge(negativeLiningTerms,1), ...
    termGlobalEdge(negativeLiningTerms,2),base(negativeLiningTerms), ...
    'VariableNames',{'elementID','node1Global','node2Global','unitKContribution'});
snapshots=table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'absoluteDay','minimumLiningTransmission_m_s','negativeLiningEdges','materiallyNegativeLiningEdges'});
if allLiningElementTermsNonnegative
    liningStatement='All individual lining element-edge contributions are nonnegative: any nonnegative element K combination has nonnegative lining transmission. This proves the coefficient sign for the entire history, conditional on the nonnegative constitutive K law.';
else
    liningStatement='Some individual lining element-edge contributions are negative: unit-K aggregation cannot prove history-wide heterogeneous-K lining positivity. Only final and explicitly saved snapshot states are evaluated; unstored accepted states remain unproved.';
    if isfield(out,'snapshots')
        for s=1:numel(out.snapshots)
            snap=out.snapshots(s);nu=2*ns;
            delta=reshape(snap.U(1:nu)-out.mechanics.referenceU(1:nu),2,[])';
            epsElement=Gx*delta(:,1)+Gz*delta(:,2);
            law=thm_bulk_constitutive(P.P*snap.T(nodes),P.P*snap.h, ...
                epsElement(P.element),out.parameters,[],false);
            sampleK=P.elementAverage*law.K;
            sampleTrans=accumarray(map,base.*repmat(sampleK,3,1),[nedge 1]);
            sampleRound=64*eps*accumarray(map,abs(base).*repmat(sampleK,3,1),[nedge 1]);
            snapshots(end+1,:)={snap.timeDay,min(sampleTrans(lining)), ...
                sum(lining&sampleTrans<0),sum(lining&sampleTrans<-sampleRound)}; %#ok<AGROW>
        end
    end
end
positive=Ke>0;positivePairs=[local(positive,[1 2]);local(positive,[2 3]);local(positive,[3 1])];
adjacency=sparse(positivePairs(:,1),positivePairs(:,2),ones(size(positivePairs,1),1),ns,ns);
component=conncomp(graph(spones(adjacency+adjacency')));
far=[];
if out.parameters.hydroFar
    q=xy(nodes,:);far=find(abs(q(:,1)-min(xy(:,1)))<1e-8 ...
        |abs(q(:,1)-max(xy(:,1)))<1e-8|abs(q(:,2)-max(xy(:,2)))<1e-8);
end
anchored=all(ismember(unique(component),unique(component(far))));
d=struct('year',year,'minimumElementK_m_s',min(Ke),'maximumElementK_m_s',max(Ke), ...
    'zeroElementKCount',sum(Ke==0),'geometricNegativeEdgeCount',sum(geometricNegative), ...
    'finalNegativeEdgeCount',sum(finalNegative),'minimumFinalSignedTransmission_m_s',min(trans), ...
    'acceptedMaximumNegativeEdgeCount',max(history(:,27)), ...
    'acceptedMinimumSignedTransmission_m_s',min(history(:,26)), ...
    'negativeLiningGeometricEdgeCount',sum(geometricNegative&lining), ...
    'materiallyNegativeLiningFinalEdgeCount',sum(lining&trans<-edgeArithmetic), ...
    'allLiningSingleElementContributionsNonnegative',allLiningElementTermsNonnegative, ...
    'historyWideLiningNonnegativityProvenForAnyNonnegativeKe',allLiningElementTermsNonnegative, ...
    'negativeLiningSingleElementContributions',liningNegativeContributions, ...
    'savedSnapshotLiningAudit',snapshots,'liningPositivityStatement',liningStatement, ...
    'diffusiveEnergy',elementEnergy,'factorizationInfinityError',factorError, ...
    'constantPatchInfinityError',constantError,'subfaceOperatorInfinityError',operatorError, ...
    'subfaceDivergenceInfinityError',divergenceError,'positiveKConnectedComponents',numel(unique(component)), ...
    'everyPositiveKComponentDirichletAnchored',anchored,'MMatrixGuaranteed',~any(finalNegative), ...
    'discreteMaximumPrincipleGuaranteed',false, ...
    'negativeGeometricEdges',edgeTable(globalEdge(geometricNegative,:),geometric(geometricNegative),trans(geometricNegative),lining(geometricNegative),xy), ...
    'negativeFinalDiffusionEdges',edgeTable(globalEdge(finalNegative,:),geometric(finalNegative),trans(finalNegative),lining(finalNegative),xy), ...
    'acceptedSignedEdgeHistory',table(history(:,1),history(:,26),history(:,27), ...
        'VariableNames',{'absoluteDay','minSignedFEMTransmission_m_s','negativeAggregateFEMEdges'}), ...
    'coercivityStatement','Nonnegative element energies establish positive semidefiniteness; strict coercivity requires anchoring every connected positive-K component. Signed cotangent coefficients do not establish a maximum principle.', ...
    'stateCoverage','Independent element/energy/flux audit is at year end. Every accepted step retains signed-edge count/min; per-step element K and all edge locations are not stored.', ...
    'negativeCoefficientInterpretation','Equivalent aggregate P1 diffusion/cross-diffusion weights, not negative physical permeability; no zero clipping was performed.');
end

function t=edgeTable(edge,geometry,transmission,lining,xy)
first=xy(edge(:,1),:);second=xy(edge(:,2),:);
t=table(edge(:,1),edge(:,2),first(:,1),first(:,2),second(:,1),second(:,2),geometry,transmission,lining, ...
    'VariableNames',{'node1Global','node2Global','x1_m','z1_m','x2_m','z2_m', ...
    'geometricWeight','signedTransmission_m_s','liningSupportEdge'});
end

function checks=inspectAggregate(A,leak,year,o)
checks=emptyChecks(); a=A.annual(A.annual(:,1)==year,:); m=A.monthly(A.monthly(:,1)==year,:);
j=A.byJoint(A.byJoint(:,1)==year,:); d=A.daily(A.daily(:,1)==year,:);
checks(end+1)=item('aggregate_annual_closure',size(a,1)==1 ...
    &&near(a(1,2:4),[leak.totalVolume_m3_per_m,sum(leak.volumeByJoint_m3_per_m),leak.bodyVolume_m3_per_m],o.volumeAbsTolerance_m2,o.relativeTolerance), ...
    'Registered annual leakage equals the original year output.');
checks(end+1)=item('aggregate_monthly_closure',size(m,1)==12&&isequal(m(:,2),(1:12)') ...
    &&near(m(:,3),leak.monthlyTotalVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(m(:,4:11),leak.monthlyVolumeByJoint_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance), ...
    'Registered months are complete and match every joint.');
checks(end+1)=item('aggregate_joint_closure',size(j,1)==8&&isequal(j(:,2),(1:8)') ...
    &&near(j(:,3),leak.jointAnglesRad(:),1e-12,o.relativeTolerance) ...
    &&near(j(:,4),leak.volumeByJoint_m3_per_m(:),o.volumeAbsTolerance_m2,o.relativeTolerance), ...
    'Registered eight-joint annual values and positions match the year output.');
checks(end+1)=item('aggregate_daily_closure',size(d,1)==365&&isequal(d(:,2),(1:365)') ...
    &&near(d(:,3),leak.dailyBinEndDay,o.timeTolerance_day,0) ...
    &&near(d(:,4),leak.dailyTotalVolume_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance) ...
    &&near(d(:,5:12),leak.dailyVolumeByJoint_m3_per_m,o.volumeAbsTolerance_m2,o.relativeTolerance), ...
    'Registered daily absolute times and all eight joint volumes match the year output.');
end

function bins=acceptedBins(t,q,edges)
% Exact integration of piecewise-constant backward-Euler accepted rates.
demand(all(diff(t)>0)&&edges(1)>=t(1)-1e-6&&edges(end)<=t(end)+1e-6, ...
    'Integration edges extend outside the complete accepted history.');
edges(1)=t(1); edges(end)=t(end);
cumulative=[zeros(1,size(q,2));cumsum(q.*(diff(t)*86400),1)];
bins=diff(interp1(t,cumulative,edges,'linear'),1,1);
demand(all(isfinite(bins(:))),'Daily/monthly integration produced nonfinite volumes.');
end

function [soilV,thermalV,orientation]=nodalVolumes(out)
xy=out.pAll; N=size(xy,1); t=out.mesh.tSoil;
area=triAreas(xy,t); soilV=accumarray(t(:),repmat(area/3,3,1),[N 1]);
thermalV=soilV;
orientation=struct('positiveElements',0,'negativeElements',0);
for name={'qTunnelGlobal','qInvertGlobal'}
    if ~isfield(out.mesh,name{1}),continue;end
    q=out.mesh.(name{1});
    for e=1:size(q,1)
        nodes=q(e,:); points=xy(nodes,:); v=zeros(4,1); gauss=[-1 1]/sqrt(3);signs=zeros(4,1);slot=0;
        for xi=gauss
            for eta=gauss
                shape=.25*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
                gradient=.25*[-(1-eta),-(1-xi);(1-eta),-(1+xi);(1+eta),(1+xi);-(1+eta),(1-xi)];
                jacobian=gradient'*points;
                determinant=det(jacobian);slot=slot+1;signs(slot)=sign(determinant);
                demand(isfinite(determinant)&&abs(determinant)>64*eps*max(norm(jacobian,'fro')^2,realmin), ...
                    'Q4 thermal element has a zero or numerically singular Gauss Jacobian.');
                v=v+shape*abs(determinant);
            end
        end
        demand(all(signs==signs(1)),'Folded Q4 thermal element changes Jacobian sign within the element.');
        orientation.positiveElements=orientation.positiveElements+(signs(1)>0);
        orientation.negativeElements=orientation.negativeElements+(signs(1)<0);
        thermalV(nodes)=thermalV(nodes)+v;
    end
end
if isfield(out.mesh,'tInvertTriGlobal')&&~isempty(out.mesh.tInvertTriGlobal)
    t=out.mesh.tInvertTriGlobal; area=triAreas(xy,t);
    thermalV=thermalV+accumarray(t(:),repmat(area/3,3,1),[N 1]);
end
soilV=soilV(out.soilNodes);
demand(all(soilV>0)&&all(thermalV>0),'Zero/negative nodal CVFE/thermal volume.');
end

function area=triAreas(xy,t)
a=xy(t(:,2),:)-xy(t(:,1),:); b=xy(t(:,3),:)-xy(t(:,1),:);
area=abs(a(:,1).*b(:,2)-a(:,2).*b(:,1))/2;
demand(all(isfinite(area))&&all(area>0),'Degenerate soil/invert triangle.');
end

function [V,baseV]=independentSlitVolume(m,p)
g=m.openingGeometry; need(g,{'nodalRadius_m','nodalGap_m'},'exact slit geometry');
r=g.nodalRadius_m; gap=g.nodalGap_m;
demand(finite(r)&&finite(gap)&&size(r,1)==8&&isequal(size(r),size(gap)), ...
    'Exact slit geometry must be finite, joint-by-radial-node data.');
demand(near(gap',m.contact.jointGap_m,1e-12,1e-10),'Geometric gaps differ from the true mechanical nodal jump.');
jointLaw=joint_microchannel_damage(p.leakage);
f=expand(p.thm.openingTransmissionFraction,8); c=expand(jointLaw.jointConductance,8); V=zeros(8,1);
for j=1:8
    [radius,order]=sort(r(j,:)); length=diff(radius); u=gap(j,order);
    demand(all(length>0),'Distinct actual mechanical radial coordinates are required.');
    if c(j)==0,continue;end
    for k=1:numel(length)
        a=u(k); b=u(k+1);
        if a>=0&&b>=0,v=length(k)*(a+b)/2;
        elseif a<=0&&b<=0,v=0;
        else,positive=max(a,b);negative=min(a,b);v=length(k)*positive^2/(2*(positive-negative));end
        V(j)=V(j)+f(j)*v;
    end
end
[~,~,inventory]=joint_channel_inventory(zeros(8,jointLaw.jointFreeze.quadraturePoints),jointLaw.jointFreeze,jointLaw.jointStorage, ...
    c,p.liningThickness,p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
baseV=inventory.waterVolume_m2;
end

function yes=eightBlocks(M)
need(M,{'q','p','nodeSegmentID','elementSegmentID','connectedComponentID','jointNodePairs'},'mechanical mesh');
nodes=M.nodeSegmentID(M.q);
yes=numel(unique(M.connectedComponentID))==8&&isequal(unique(M.nodeSegmentID),(1:8)') ...
    &&all(nodes==M.elementSegmentID,'all')&&size(M.jointNodePairs,3)==8;
for j=1:8
    pairs=M.jointNodePairs(:,:,j);
    yes=yes&&all(pairs(:,1)~=pairs(:,2))&&near(M.p(pairs(:,1),:),M.p(pairs(:,2),:),1e-10,0) ...
        &&all(M.nodeSegmentID(pairs(:,1))~=M.nodeSegmentID(pairs(:,2)));
end
end

function y=year20Summary(out,potential,exactV,baseV)
leak=out.leakage; h=out.history;
y=struct('year20AnnualLeakage_L_per_m',1000*leak.totalVolume_m3_per_m, ...
    'year20ByJoint_L_per_m',1000*leak.volumeByJoint_m3_per_m(:), ...
    'year20ByMonth_L_per_m',1000*leak.monthlyTotalVolume_m3_per_m(:), ...
    'dailyActual_L_per_day_per_m',1000*leak.dailyTotalVolume_m3_per_m(:), ...
    'monthlyActualByJoint_L_per_m',1000*leak.monthlyVolumeByJoint_m3_per_m, ...
    'sameStateFarHead_L_per_m',1000*sum(potential.annual), ...
    'sameStateFarHeadByJoint_L_per_m',1000*potential.annual(:), ...
    'sameStateFarHeadByMonth_L_per_m',1000*sum(potential.monthly,2), ...
    'sameStateFarHeadMonthlyByJoint_L_per_m',1000*potential.monthly, ...
    'sameStateFarHeadDaily_L_per_day_per_m',1000*sum(potential.daily,2), ...
    'supplyDifference_L_per_m',1000*(sum(potential.annual)-leak.totalVolume_m3_per_m), ...
    'supplyLossFraction',divideDefined(sum(potential.annual)-leak.totalVolume_m3_per_m,sum(potential.annual)), ...
    'differenceDenominator','QsameStateFarHead; negative difference is allowed for stored/thermal-pressure discharge.', ...
    'comparisonDefinition','Original accepted total Ceff(T,u) times quadrature-positive FAR head after entry; instantaneous T and opening are identical to the full solve.', ...
    'finalExactSlitGeometricVolume_m2',exactV,'baseChannelWaterVolume_m2',baseV);
open=out.mechanics.jointHistory(:,58:65);total=out.mechanics.jointHistory(:,66:73);
share=zeros(size(open));positive=total>0;share(positive)=open(positive)./total(positive);
y.maxOpenSlitConductance_m_s=max(open,[],'all');
y.openSlitActiveDuration_day=sum(h(any(open>0,2),2));
y.maxOpenConductanceShare=max(share,[],'all');
y.maxCapillaryConductance_m_s=max(out.mechanics.jointHistory(:,50:57),[],'all');
y.openSlitInterpretation='An implemented/independently tested branch is distinct from activation in this actual history; no gasket damage/failure is inferred.';
actual=sum(leak.jointRate_m2_per_s,2); simple=sum(potential.rate_m2_s,2);
ratio=NaN(size(actual)); positive=simple>0; ratio(positive)=actual(positive)./simple(positive);
y.acceptedHydraulics=table(h(:,1)-19*365,h(:,2),actual,simple,ratio, ...
    'VariableNames',{'dayOfYear','dt_day','actualRate_m2_s','sameStateFarHeadRate_m2_s','actualOverSameStateFarHead'});
y.jointAngles_rad=leak.jointAnglesRad(:);
end

function m=metricTemplate()
names={'year','acceptedSteps','rejectedSteps','unknowns','duration_day','totalLeak_m3_per_m', ...
    'jointLeak_m3_per_m','bodyLeak_m3_per_m','sameStateFarHeadLeak_m3_per_m', ...
    'supplyDifference_m3_per_m','supplyLossFraction','waterSignedAnnualDefect_m2', ...
    'energySignedAnnualDefect_J_per_m','waterAbsoluteAnnualDefect_m2','energyAbsoluteAnnualDefect_J_per_m', ...
    'maxGlobalWaterDefect_m2','maxGlobalEnergyDefect_J_per_m','maxLocalWaterDefect_m2', ...
    'maxLocalEnergyDefect_J_per_m','globalWaterTolerance_m2','localWaterTolerance_m2', ...
    'globalEnergyTolerance_J_per_m','localEnergyTolerance_J_per_m','maxScaledResidual', ...
    'maxMechanicalResidual_N_per_m','mechanicalTolerance_N_per_m','maxIncrementalVolumetricStrain', ...
    'maxIncrementalPrincipalStrain','maxIncrementalTensorShearStrain','maxIncrementalGradientNorm', ...
    'maxAbsolutePrincipalStrain','maxAbsoluteGradientNorm','installedPrincipalStrain', ...
    'installedGradientNorm','maxOmittedSlitVolumeFraction','minFinalPorosity','maxFinalPorosity', ...
    'maxOpenSlitConductance_m_s','maxCapillaryConductance_m_s','openSlitActiveDuration_day', ...
    'maxOpenConductanceShare','positiveQ4Elements','negativeQ4Elements'};
m=cell2struct(num2cell(NaN(1,numel(names))),names,2);
end

function s=parameterSummary(p,out)
[jointLaw,jointDamage]=joint_microchannel_damage(p.leakage);
s=struct('status','Illustrative parameters; no field or joint-test calibration supplied.', ...
    'soilKsat_m_s',p.Ksat,'initialPorosity',p.porosity,'soilElasticity_Pa',p.mechanics.Esoil, ...
    'soilPoissonRatio',p.mechanics.nuSoil,'soilBiot',p.mechanics.biot, ...
    'concreteElasticity_Pa',p.mechanics.Econc,'liningK_m_s',p.liningK, ...
    'jointConductance_m_s',p.leakage.jointConductance,'jointEntryHead_m',p.leakage.jointEntryHead, ...
    'contactPenalty_Pa_per_m',p.thm.contactNormalStiffness_Pa_per_m, ...
    'boltStiffness_N_per_m2',p.thm.boltNormalStiffness_N_per_m2, ...
    'boltPreload_N_per_m',p.thm.boltPreload_N_per_m, ...
    'openingTransmissionFraction',p.thm.openingTransmissionFraction, ...
    'waterTableDepth_m',p.waterTableDepth,'forcing',p.forcing,'initial',p.initial, ...
    'topAirRange_C',[p.topMin p.topMax],'tunnelAirRange_C',[p.tunnelMin p.tunnelMax], ...
    'topHTC_W_m2_K',p.topHTC,'tunnelHTC_W_m2_K',p.tunnelHTC, ...
    'geothermalFlux_W_m2',p.geothermalFlux,'soilViscosity',p.soilViscosity, ...
    'jointFreeze',p.leakage.jointFreeze,'jointStorage',p.leakage.jointStorage, ...
    'thm',p.thm,'solver',p.solver,'thermalNodes',size(out.pAll,1), ...
    'soilNodes',numel(out.soilNodes),'soilTriangles',size(out.mesh.tSoil,1), ...
    'ringMechanicalNodes',size(out.mechanics.mesh.p,1),'coupledUnknowns',numel(out.finalX));
if isfield(p,'meshFile'),s.meshFile=p.meshFile;end
s.resolvedJointDamage=jointDamage;
s.resolvedJointConductance_m_s=jointLaw.jointConductance;
s.resolvedJointFreeze=jointLaw.jointFreeze;
end

function figs=plotReport(r,o)
y=r.yearSummary; q=r.year20; figs=struct();
figs.annual=figure('Visible',o.visible,'Color','w','Position',[100 100 860 390]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile; plot(ax,y.year,1000*y.totalLeak_m3_per_m,'o-','LineWidth',1.3,'MarkerSize',4);hold(ax,'on');
plot(ax,y.year,1000*y.sameStateFarHeadLeak_m3_per_m,'s--','LineWidth',1.2,'MarkerSize',4);
xlabel(ax,'Simulated year');ylabel(ax,'Annual leakage (L m^{-1})');title(ax,'(a) Complete coupled history');
legend(ax,{'Coupled THM','Same-state C_{eff}, far head'},'Location','best','Box','off');style(ax);
ax=nexttile; plot(ax,y.year,100*y.supplyLossFraction,'o-','LineWidth',1.3,'MarkerSize',4);
yline(ax,0,':');xlabel(ax,'Simulated year');ylabel(ax,'Supply difference / Q_{far} (%)');
title(ax,'(b) Supply restriction at the same state');style(ax);
figs.seasonal=figure('Visible',o.visible,'Color','w','Position',[100 100 860 390]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile;plot(ax,1:365,q.dailyActual_L_per_day_per_m,'-','LineWidth',1.2);hold(ax,'on');
plot(ax,1:365,q.sameStateFarHeadDaily_L_per_day_per_m,'--','LineWidth',1.2);
xlabel(ax,'Day in year 20');ylabel(ax,'Daily leakage (L d^{-1} m^{-1})');title(ax,'(a) Conservative daily bins');style(ax);
legend(ax,{'Coupled THM','Same-state C_{eff}, far head'},'Location','best','Box','off');
ax=nexttile;bar(ax,(1:12)',[q.year20ByMonth_L_per_m,q.sameStateFarHeadByMonth_L_per_m],'grouped');
xlabel(ax,'Month in year 20');ylabel(ax,'Monthly leakage (L m^{-1})');xticks(ax,1:12);
title(ax,'(b) Equal month-bin volumes');style(ax);
figs.joints=figure('Visible',o.visible,'Color','w','Position',[100 100 860 390]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile;bar(ax,(1:8)',[q.year20ByJoint_L_per_m,q.sameStateFarHeadByJoint_L_per_m],'grouped');
xlabel(ax,'Longitudinal joint');ylabel(ax,'Year-20 leakage (L m^{-1})');xticks(ax,1:8);
specifiedJoints=find(expand(r.parameterSummary.jointConductance_m_s,8)>0);
if isempty(specifiedJoints), jointCaption='All eight joints prescribed sealed';
else, jointCaption=['Prescribed defect joint(s): ' strtrim(sprintf('%d ',specifiedJoints))]; end
title(ax,{'(a) Computed joint leakage',jointCaption});style(ax);
legend(ax,{'Coupled THM','Same-state C_{eff}, far head'}, ...
    'Location','northwest','Box','off','FontSize',9);
ax=nexttile;imagesc(ax,1:12,1:8,q.monthlyActualByJoint_L_per_m');set(ax,'YDir','normal');
xlabel(ax,'Month in year 20');ylabel(ax,'Longitudinal joint');xticks(ax,1:12);yticks(ax,1:8);
title(ax,'(b) Computed joint-month volumes');style(ax);grid(ax,'off');
c=colorbar(ax);c.Label.String='Leakage (L m^{-1})';colormap(ax,parula);
end

function style(ax)
set(ax,'FontName','Arial','FontSize',10,'LineWidth',.8,'Box','on');grid(ax,'on');ax.GridAlpha=.15;
end

function exportReport(report,figs,folder)
writetable(report.checks,fullfile(folder,'thm_20year_checks.csv'));
writetable(report.yearSummary,fullfile(folder,'thm_20year_year_summary.csv'));
for k=1:numel(report.hydraulicAudit)
    audit=report.hydraulicAudit(k);prefix=sprintf('thm_year%02d_',audit.year);
    writetable(audit.negativeGeometricEdges,fullfile(folder,[prefix 'negative_geometric_edges.csv']));
    writetable(audit.negativeFinalDiffusionEdges,fullfile(folder,[prefix 'negative_final_FEM_edges.csv']));
    writetable(audit.acceptedSignedEdgeHistory,fullfile(folder,[prefix 'signed_FEM_edge_history.csv']));
    writetable(audit.negativeLiningSingleElementContributions,fullfile(folder,[prefix 'negative_lining_element_contributions.csv']));
    writetable(audit.savedSnapshotLiningAudit,fullfile(folder,[prefix 'saved_snapshot_lining_audit.csv']));
end
if report.passed
    q=report.year20;
    writetable(table((1:12)',q.year20ByMonth_L_per_m,q.sameStateFarHeadByMonth_L_per_m, ...
        'VariableNames',{'month','actual_L_m','sameStateFarHead_L_m'}),fullfile(folder,'thm_year20_monthly.csv'));
    writetable(table((1:8)',q.jointAngles_rad,q.year20ByJoint_L_per_m,q.sameStateFarHeadByJoint_L_per_m, ...
        'VariableNames',{'jointID','angle_rad','actual_L_m','sameStateFarHead_L_m'}),fullfile(folder,'thm_year20_joint.csv'));
    writetable(q.acceptedHydraulics,fullfile(folder,'thm_year20_accepted_hydraulics.csv'));
    writetable(array2table(q.monthlyActualByJoint_L_per_m,'VariableNames', ...
        arrayfun(@(j)sprintf('joint%d_L_m',j),1:8,'UniformOutput',false)),fullfile(folder,'thm_year20_monthly_by_joint.csv'));
end
for key=fieldnames(figs)'
    handle=figs.(key{1});
    if isgraphics(handle,'figure')
        exportgraphics(handle,fullfile(folder,['thm_20year_' key{1} '.png']),'Resolution',300);
        exportgraphics(handle,fullfile(folder,['thm_20year_' key{1} '.pdf']),'ContentType','vector');
    end
end
fid=fopen(fullfile(folder,'thm_20year_validation.txt'),'w','n','UTF-8');
demand(fid>=0,'Cannot open validation text report.');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'THM consecutive-history numerical audit\npassed=%d; complete=%g/20; schema=%g; engine=%s\n', ...
    report.passed,report.completedYears,report.schemaVersion,report.engineVersion);
fprintf(fid,'Physical field validation: none. Source SHA-256 is computed at validation time; the annual driver did not record an immutable runtime source hash.\n');
if report.passed
    fprintf(fid,'Year20 actual=%.12g L/m; same-state far-head=%.12g L/m; difference/Qfar=%.9g%%.\n', ...
        report.year20.year20AnnualLeakage_L_per_m,report.year20.sameStateFarHead_L_per_m,100*report.year20.supplyLossFraction);
end
for k=1:numel(report.issues),fprintf(fid,'FAIL %s\n',report.issues{k});end
for k=1:numel(report.warnings),fprintf(fid,'NOTE %s\n',report.warnings{k});end
for k=1:numel(report.limitations),fprintf(fid,'LIMIT %s\n',report.limitations{k});end
for k=1:numel(report.sources),fprintf(fid,'SOURCE %s | %s\n',report.sources(k).doi,report.sources(k).url);end
fprintf(fid,'Configuration JSON:\n%s\n',jsonencode(report.configuration));
fprintf(fid,'Physical parameter summary JSON:\n%s\n',jsonencode(report.parameterSummary));
fprintf(fid,'Inspected source SHA-256 JSON:\n%s\n',jsonencode(report.sourceHash));
fprintf(fid,'Optional campaign source-at-start snapshot SHA-256 JSON:\n%s\n',jsonencode(report.sourceAtStartHash));
save(fullfile(folder,'thm_20year_validation.mat'),'report','-v7');
end

function h=sourceHashes(root)
if nargin<1,root=fileparts(mfilename('fullpath'));end
names={'TH_2D_segment_v16.m','run_TH_segment_year20_thm.m', ...
    'joint_contact_thm_model.m','joint_supply_th_model.m','soil_poroelastic_fem.m', ...
    'solve_segment_contact_fem.m','thm_bulk_constitutive.m','frozen_soil_constitutive.m', ...
    'build_soil_cvfe_flux.m','soil_cvfe_flux.m','soil_cvfe_heat.m','build_soil_scv_quadrature.m', ...
    'joint_opening_profile_hydraulics.m','joint_opening_p1_hydraulics.m', ...
    'joint_profile_geometric_volume.m','evaluate_joint_drain.m','validate_thm_20year.m', ...
    'joint_microchannel_damage.m','thm_constitutive_signature.m', ...
    'joint_cycle_damage_defaults.m','joint_cycle_channel_law.m', ...
    'joint_freeze_thaw_memory.m','joint_cycle_severity.m', ...
    'validate_joint_cycle_memory.m','audit_joint_cycle_damage.m','audit_joint_channel_geometry.m'};
h=struct('file',{},'sha256',{},'status',{});
for k=1:numel(names)
    file=fullfile(root,names{k}); digest=''; status='missing';
    try
        if isfile(file)
            fid=fopen(file,'rb'); demand(fid>=0,'Cannot read source hash file.');
            cleanup=onCleanup(@()fclose(fid)); bytes=fread(fid,inf,'*uint8');
            md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes(:),'int8'));
            raw=typecast(md.digest(),'uint8');digest=lower(reshape(dec2hex(raw,2).',1,[]));
            clear cleanup;status='current source at validation time; not a recorded runtime fingerprint';
        end
    catch ME,status=ME.message;end
    h(end+1)=struct('file',file,'sha256',digest,'status',status); %#ok<AGROW>
end
end

function s=researchSources()
doi={'10.1016/j.tust.2025.106492','10.1016/j.compgeo.2025.107186', ...
    '10.1016/j.tust.2026.107519','10.1016/j.compgeo.2026.108448','10.1016/j.enggeo.2025.108136', ...
    '10.1016/j.compgeo.2026.108114'};
url={'https://www.sciencedirect.com/science/article/pii/S0886779825001300', ...
    'https://repository.tudelft.nl/record/uuid%3A1739ae5d-3d49-4c21-8093-678ec2d9b978', ...
    'https://www.sciencedirect.com/science/article/pii/S0886779826000775', ...
    'https://publications.polymtl.ca/80436/', ...
    'https://www.sciencedirect.com/science/article/pii/S0013795225002327', ...
    'https://publications.rwth-aachen.de/record/1035539/files/1035539.pdf'};
role={'Gasket compression/deformation and leakage require independent joint experiments; no gasket failure law imported.', ...
    'Independent-pressure triple-node cohesive THM interface; the present Robin/contact model is not that interface.', ...
    'Long-term cold-region tunnel THM with experimental rock freeze-thaw degradation; no calibrated rock or gasket mechanical damage law is imported.', ...
    'Frozen-soil TEVP/free-energy framework; primary public abstract/author record only, no full text acquired, no TEVP implemented.', ...
    'Multi-observation soil-column freeze-thaw calibration including localized ice/strain; present model has no ice lenses/hysteretic residual settlement.', ...
    'Full public paper inspected; equations 18, 20-23 and 25-26 reproduced at material-point level only. No Brazilian FE reproduction or transfer of rock fatigue parameters to gaskets. The optional local-cycle hydraulic coarsening law is a separate uncalibrated research hypothesis.'};
s=struct('doi',doi,'url',url,'relationship',role);
end

function c=limitations()
c={'Numerical audit only: no field data or laboratory parameter calibration was supplied.', ...
    'A passing history is not mesh/time-step convergence and does not establish the model accuracy for reliability analysis.', ...
    'Obtuse/heterogeneous P1 diffusion can have negative aggregate cotangent edge weights despite nonnegative element permeability and positive energy. No M-matrix/discrete-maximum-principle guarantee is claimed; weights are never clipped.', ...
    'Equal water/ice density inventory: no true 9% density expansion or segregated ice lenses.', ...
    'Elastic soil/ring/contact: no soil plasticity, TEVP creep or calibrated mechanical damage evolution. Fixed joint-channel damage is a sensitivity scenario; optional local-cycle constant-volume hydraulic coarsening is an illustrative uncalibrated research hypothesis.', ...
    'Local-cycle coarsening updates hydraulic geometry one accepted step after thaw completion; it is sequential, not a fully implicit damage/THM solve. Connected void volume is held fixed; surface-energy and damage-heat terms are omitted.', ...
    'Eight mechanical ring blocks exclude invert/slab mechanics; heat mesh stays continuous across joints.', ...
    'Soil-lining interface is a bilateral linear spring coupling; no interface separation/sliding law.', ...
    'Fixed reference geometry for transport and gravity: no ALE/deformation-gradient pullback; omitted geometry is O(strain).', ...
    'Existing-defect Robin flow only; no hydraulic fracture, calibrated gasket failure or direct distributed slit-pressure traction.', ...
    'New open slit uses zero storage: the exact <=1% geometric inventory guard is instantaneous, not a bound on annual leakage error.', ...
    'Mechanical work is omitted from the heat equation; research references do not validate these omitted terms.', ...
    'Inspected current source hashes cannot reconstruct an immutable runtime code hash absent from the stored annual schema.'};
end

function c=emptyChecks()
c=struct('name',{},'passed',{},'detail',{});
end
function c=item(name,passed,detail)
c=struct('name',name,'passed',logical(passed),'detail',detail);
end
function need(s,names,label)
demand(isstruct(s)&&isscalar(s)&&all(isfield(s,names)), ...
    sprintf('Missing required fields in %s: %s.',label,strjoin(names,', ')));
end
function demand(condition,message)
assert(isscalar(condition)&&logical(condition),'THM20Audit:Data','%s',message);
end
function yes=finite(x)
yes=isnumeric(x)&&isreal(x)&&~isempty(x)&&all(isfinite(x(:)));
end
function yes=near(a,b,absolute,relative)
yes=isnumeric(a)&&isnumeric(b)&&isreal(a)&&isreal(b)&&isequal(size(a),size(b)) ...
    &&all(isfinite(a(:)))&&all(isfinite(b(:))) ...
    &&all(abs(a(:)-b(:))<=absolute+relative*max(abs(a(:)),abs(b(:))));
end
function x=expand(x,n)
demand(finite(x)&&isvector(x)&&(isscalar(x)||numel(x)==n),'Invalid scalar/per-joint parameter.');
x=double(x(:));if isscalar(x),x=repmat(x,n,1);end
end
function value=divideDefined(numerator,denominator)
value=NaN;if denominator>0,value=numerator/denominator;end
end
function yes=isSwitch(x)
yes=(isnumeric(x)||islogical(x))&&isreal(x)&&isscalar(x)&&isfinite(x)&&ismember(x,[0 1]);
end
