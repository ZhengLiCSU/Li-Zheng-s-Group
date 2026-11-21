function out = TH_2D_segment_v16(opt)
% TH_2D_SEGMENT_V16 Conservative monolithic THM with eight contact segments.
% Positive SCV phase quadrature and physical Darcy subface heat advection.
%   out = TH_2D_segment_v16(opt) advances one tunnel realization.
%   opt supplies the mesh, soil and lining properties, initial fields and time-step controls.
%   The accepted-step records are in out.history; out.historyColumns names each column.
%   Qout_m2s is the outflow per metre of tunnel axis. The leakage summary is in out.leakage.
% T [degC], h [m LIQUID pressure head], p=rho*g*h [Pa], Z positive down.
% Painter & Karra (2014), Eqs 10-16, beta=omega=1 special case; VG-Mualem.
% Biot pore-volume feedback with an equal ice/water density approximation.
% Monolithic small-strain THM with pore-ice pressure and eight contact blocks.
% No independent moisture PDE, pressure clipping, or artificial relaxation.
% Joint conductance and hydraulic head-loss threshold require calibration.
if nargin<1, opt=struct(); end
if ischar(opt) && strcmp(opt,'selftest'), out=selftest(); return; end
if ischar(opt) && strcmp(opt,'thmselftest'), out=thmSelftest(); return; end
p=defaults(); p=merge(p,opt); validateParameters(p);
if isempty(p.mesh)
    if ~isfile(p.meshFile)
        error('Run the supplied mesh script first, or set opt.meshFile to its MAT output.');
    end
    S=load(p.meshFile);
else
    S=p.mesh;
end
G=geometry(S,p); B=boundaries(G,S,p);
if p.leakage.jointFreeze.enabled
    G=jointCouplingPattern(G,B);
    G.jointThermalPaths=B.jointThermalPaths;
    Q=B.jointQuadrature;
    nj=Q.nj; nq=B.jointThermalPaths.nq;
    Jmap=sparse((1:numel(Q.weight))',Q.jointID,1,numel(Q.weight),nj);
    Wj=Jmap'*spdiags(Q.weight,0,numel(Q.weight),numel(Q.weight))*Q.phi;
    [~,~,inventory]=joint_channel_inventory(zeros(nj,nq),G.jointLaw.jointFreeze, ...
        G.jointLaw.jointStorage,G.jointLaw.jointConductance,p.liningThickness, ...
        p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
    G.jointChannelGeometry=audit_joint_channel_geometry(inventory, ...
        G.jointLaw.jointWidth_m,p.liningThickness);
    G.channelWaterNodeMass=Wj'*inventory.waterVolume_m2;
    mappedVolume=B.jointThermalPaths.interpolation'* ...
        repelem(inventory.waterVolume_m2/nq,nq);
    assert(all(mappedVolume<=G.vconc*(1+1e-10)), ...
        'TH:channelGeometry','Connected channel volume exceeds the concrete volume it replaces. Reconsider conductance/radius/tortuosity.');
else
    G.channelWaterNodeMass=zeros(G.ns,1);
end
if p.thm.enabled, G=thmGeometry(G,S,p); end
[X,t0,initInfo,G]=initialFields(G,B,p);
jointCycleHistory=zeros(0,49);jointCycleIceHistory=[];
if G.jointCycleSettings.enabled && isempty(G.jointCycleMemory)
    initialCycleState=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
    G.jointCycleMemory=joint_freeze_thaw_memory([], ...
        1-initialCycleState.jointInventory.liquidVolumeFraction,t0/p.day, ...
        G.fixedJointLaw.jointConductance(:)>0,G.jointCycleSettings);
end
initInfo.jointCycleInitialMemory=G.jointCycleMemory;
if p.verbose
    fprintf('TH v16: %d thermal + %d liquid-head DOFs; SCV phase quadrature and physical subface fluxes.\n',G.N,G.ns);
    fprintf('Equal-density, small-strain THM=%d; %d mechanical DOFs.\n',p.thm.enabled,G.nm);
end
start=X; mechanicalHistory=zeros(0,6); jointMechanicalHistory=zeros(0,73);
strainHistory=zeros(0,7);
hist=[]; cycles=[]; snaps=struct('timeDay',{},'T',{},'h',{},'thetaW',{},'thetaI',{},'U',{});
maxAcceptedLocalWaterDefect=0; maxAcceptedLocalEnergyDefect=0;
leakHistSegment=zeros(0,p.leakage.segmentCount);
leakHistJoint=zeros(0,numel(p.leakage.jointAnglesRad));
leakHistBody=zeros(0,1);
jointFreezeHist=zeros(0,4*numel(p.leakage.jointAnglesRad));
jointHydraulicHist=zeros(0,3*numel(p.leakage.jointAnglesRad)+3);
yearRef=X; cycleConverged=false; t=t0; dt=p.dtInitial; naccepted=0;
priorStart=[]; priorAcceptedDt=NaN; tangentCache=[];
nrejected=initInfo.priorRejectedSteps;
endTime=t0+p.durationDays*p.day; nextReport=t0+p.outputEveryDays*p.day;
nextYear=t0+p.yearDays*p.day; productionEnd=inf;
if ~strcmp(p.initial.mode,'spinup'), productionEnd=endTime; end
if strcmp(p.initial.mode,'checkpoint') && isfield(p.initial,'resumeHistory') && p.initial.resumeHistory
    saved=load(p.initial.checkpointFile);
    assert(isfield(saved.checkpoint,'runInventory'),'THM:ResumeInventory', ...
        'Checkpoint lacks complete accepted-history inventory.');
    run=saved.checkpoint.runInventory;
    if G.jointCycleSettings.enabled
        assert(isfield(run,'jointCycleHistory'),'JointCycle:Resume', ...
            'Cycle-damage recovery requires its complete accepted history.');
        jointCycleHistory=run.jointCycleHistory;
        jointCycleIceHistory=run.jointCycleIceHistory;
        initInfo.jointCycleInitialMemory=run.startCycleMemory;
    end
    assert(isnumeric(run.startX)&&numel(run.startX)==numel(X)&&all(isfinite(run.startX(:))) ...
        &&isfinite(run.startTime)&&isfinite(run.endTime)&&run.startTime<t ...
        &&isfinite(run.nextReport)&&run.nextReport>t ...
        &&isfinite(run.nextYear)&&run.nextYear>t ...
        &&numel(run.yearRef)==numel(X)&&all(isfinite(run.yearRef(:))) ...
        &&isfinite(run.nextDt)&&run.nextDt>0 ...
        &&(isempty(run.priorStart)||(numel(run.priorStart)==numel(X)&&all(isfinite(run.priorStart(:))))), ...
        'THM:ResumeInventory','Invalid restored full-state/time-step inventory.');
    assert(run.endTime>=t && run.startTime<t ...
        &&abs((run.endTime-run.startTime)-p.durationDays*p.day)<1e-6, ...
        'THM:ResumeInventory','Restored histories must match the original run duration and may not pass its end.');
    t0=run.startTime; start=run.startX; endTime=run.endTime; productionEnd=endTime;
    nextReport=run.nextReport; nextYear=run.nextYear; yearRef=run.yearRef;
    dt=run.nextDt; priorStart=run.priorStart; priorAcceptedDt=run.priorAcceptedDt;
    maxAcceptedLocalWaterDefect=run.maxLocalWaterDefect;
    maxAcceptedLocalEnergyDefect=run.maxLocalEnergyDefect;
    hist=saved.hist; cycles=saved.cycles; snaps=saved.snaps;
    leakHistSegment=saved.leakHistSegment; leakHistJoint=saved.leakHistJoint;
    leakHistBody=saved.leakHistBody; jointFreezeHist=saved.jointFreezeHist;
    jointHydraulicHist=saved.jointHydraulicHist;
    mechanicalHistory=saved.mechanicalHistory; jointMechanicalHistory=saved.jointMechanicalHistory;
    strainHistory=saved.strainHistory; naccepted=size(hist,1);
    if G.jointCycleSettings.enabled
        assert(isequal(size(jointCycleHistory),[naccepted 49]) ...
            &&all(isfinite(jointCycleHistory(:))) ...
            &&isequal(jointCycleHistory(:,1),hist(:,1)), ...
            'JointCycle:Resume','Cycle history must cover every accepted step.');
    end
    assert(~isempty(hist)&&all(isfinite(hist(:)))&&all(hist(:,2)>0) ...
        &&all(diff(hist(:,1))>0) ...
        &&max(abs(diff([t0/p.day;hist(:,1)])-hist(:,2)))<1e-6, ...
        'THM:ResumeInventory','Restored history has duplicated, missing, or nonfinite accepted steps.');
    arrays={leakHistSegment,leakHistJoint,leakHistBody,jointFreezeHist,jointHydraulicHist};
    if p.thm.enabled, arrays=[arrays,{mechanicalHistory,jointMechanicalHistory,strainHistory}]; end
    assert(all(cellfun(@(v)size(v,1)==naccepted&&all(isfinite(v(:))),arrays)), ...
        'THM:ResumeInventory','Restored diagnostic histories do not cover every accepted step.');
    if p.thm.enabled
        assert(max(abs(mechanicalHistory(:,1)-hist(:,1)))<1e-8 ...
            &&max(abs(jointMechanicalHistory(:,1)-hist(:,1)))<1e-8 ...
            &&max(abs(strainHistory(:,1)-hist(:,1)))<1e-8, ...
            'THM:ResumeInventory','Mechanical/strain history dates differ from accepted steps.');
    end
    assert(abs(hist(end,1)*p.day-t)<1e-6 && abs(sum(hist(:,2))*p.day-(t-t0))<1e-4, ...
        'THM:ResumeInventory','Restored accepted history must end exactly at the checkpoint.');
    % At an annual boundary the monthly checkpoint precedes the cycle
    % diagnostic and results MAT. Recover that interrupted finalization
    % without taking a new physical step or duplicating the annual record.
    if abs((t-t0)-p.yearDays*p.day)<1e-6 ...
            &&(isempty(cycles)||cycles(end,1)<1-1e-9)
        last=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
        first=state(start(1:G.N),start(G.N+(1:G.ns)),G,p,start(G.uIndex));
        cycles(end+1,:)=[1,max(abs(X(1:G.N)-start(1:G.N))), ...
            max(abs(X(G.N+(1:G.ns))-start(G.N+(1:G.ns)))), ...
            max(abs(last.total-first.total)),max(abs(last.thetaI-first.thetaI))];
    end
end
% This cache is local to one geometry/installed-reference/parameter set.
% It is reconstructed after restart and only advanced on accepted steps.
acceptedState=[];
while t < min(endTime,productionEnd)-1e-7
    proposedDt=dt;
    [calendarCap,calendarBoundary]=calendarStepLimit(t,p);
    target=min([t+dt,t+calendarCap,endTime,productionEnd,nextReport,nextYear,calendarBoundary]); d=target-t;
    seed=[];
    if p.solver.predictor && ~isempty(priorStart) && priorAcceptedDt>0
        extrapolation=min(d/priorAcceptedDt,4);
        seed=X+extrapolation*(X-priorStart);
        if any(~isfinite(seed)) || any(seed(1:G.N)<=-100 | seed(1:G.N)>=100)
            seed=[];
        end
    end
    stepGeometry=G;stepCache=tangentCache;stepOldState=acceptedState;
    if G.jointCycleSettings.enabled
        % Apply only damage queued by an earlier ACCEPTED thaw completion.
        % The old energy uses its old radius distribution, even with state
        % caching disabled. Rejections commit neither geometry nor memory.
        if isempty(stepOldState)
            stepOldState=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
        end
        if ~isequal(G.jointCycleAppliedSeverity,G.jointCycleMemory.targetSeverity)
            stepGeometry=applyJointCycleGeometry(G,G.jointCycleMemory.targetSeverity);
            stepCache=[];
        end
    end
    [Xtrial,info,newCache,trialState,oldState]=step(X,t,d,stepGeometry,B,p,seed,stepCache,stepOldState);
    if ~info.ok && ~isempty(seed)
        % A predictor cannot turn an otherwise solvable step into a rejection.
        [Xtrial,info,newCache,trialState,oldState]=step(X,t,d,stepGeometry,B,p,[],[],oldState);
    end
    if ~info.ok
        tangentCache=[];
        nrejected=nrejected+1; dt=d/2;
        if dt<p.dtMin
            if p.saveOutput
                checkpoint=makeCheckpoint(X,t,G,p,initInfo,nrejected); %#ok<NASGU>
                checkpoint.runInventory=runInventory();
                save([p.outputPrefix '_FAILED_checkpoint.mat'],'checkpoint','hist','cycles','snaps', ...
                    'leakHistSegment','leakHistJoint','leakHistBody','jointFreezeHist','jointHydraulicHist', ...
                    'mechanicalHistory','jointMechanicalHistory','strainHistory','info');
            end
            error('TH:nonconvergence','Step at day %.5f failed: %s. Last accepted state retained; reduce forcing/inspect parameters.',t/p.day,info.reason);
        end
        if p.verbose, fprintf('  retry day %.4f: dt -> %.5g day (%s)\n',t/p.day,dt/p.day,info.reason); end
        continue
    end
    if G.jointCycleSettings.enabled && isempty(trialState)
        trialState=state(Xtrial(1:G.N),Xtrial(G.N+(1:G.ns)), ...
            stepGeometry,p,Xtrial(G.uIndex));
    end
    if p.solver.reuseConstitutiveState || G.jointCycleSettings.enabled
        [~,D]=residual(Xtrial,X,target,d,stepGeometry,B,p,false,oldState,trialState);
        acceptedState=trialState;
    else
        [~,D]=residual(Xtrial,X,target,d,stepGeometry,B,p,false);
    end
    maxAcceptedLocalWaterDefect=max(maxAcceptedLocalWaterDefect,max(abs(D.localWaterDefect)));
    maxAcceptedLocalEnergyDefect=max(maxAcceptedLocalEnergyDefect,max(abs(D.localEnergyDefect)));
    priorStart=X; priorAcceptedDt=d;
    tangentCache=newCache;
    G=stepGeometry;
    X=Xtrial; t=target; naccepted=naccepted+1;
    if G.jointCycleSettings.enabled
        G.jointCycleMemory=joint_freeze_thaw_memory(G.jointCycleMemory, ...
            1-trialState.jointInventory.liquidVolumeFraction,t/p.day, ...
            G.fixedJointLaw.jointConductance(:)>0,G.jointCycleSettings);
        jointCycleHistory(end+1,:)=[t/p.day,G.jointCycleAppliedSeverity', ...
            G.jointCycleMemory.targetSeverity',mean(G.jointCycleMemory.cycleCount,2)', ...
            mean(G.jointCycleMemory.cumulativeDose,2)', ...
            trialState.jointInventory.waterVolume_m2',G.jointLaw.jointConductance(:)']; %#ok<AGROW>
        jointCycleIceHistory(end+1,:)=reshape( ...
            (1-trialState.jointInventory.liquidVolumeFraction)',1,[]); %#ok<AGROW>
    end
    if p.thm.enabled
        mechanicalHistory(end+1,:)=[t/p.day,D.maxVolumetricStrain,D.maxJointGap_m, ...
            D.maxJointPressure_Pa,D.mechanicalResidual_Nm,D.maxOmittedSlitVolumeFraction]; %#ok<AGROW>
        jointMechanicalHistory(end+1,:)=D.jointMechanicalRecord; %#ok<AGROW>
        strainHistory(end+1,:)=D.strainRecord; %#ok<AGROW>
    end
    rec=[t/p.day,d/p.day,info.iter,D.water,D.energy,D.Qout,D.Eout,D.waterError,D.energyError,D.leak,info.norm, ...
        info.jacobianBuilds,info.linearSolves,info.residualEvaluations];
    hist=[hist;rec]; %#ok<AGROW>
    leakHistSegment=[leakHistSegment;D.leakBySegment(:)']; %#ok<AGROW>
    leakHistJoint=[leakHistJoint;D.leakByJoint(:)']; %#ok<AGROW>
    leakHistBody=[leakHistBody;D.leakBody]; %#ok<AGROW>
    jointFreezeHist=[jointFreezeHist;D.jointFreezeRecord]; %#ok<AGROW>
    jointHydraulicHist=[jointHydraulicHist;D.jointHydraulicRecord]; %#ok<AGROW>
    dt=acceptedTimeStep(proposedDt,d,info,p);
    if t>=nextReport-1e-7
        a=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
        snaps(end+1)=struct('timeDay',t/p.day,'T',X(1:G.N),'h',X(G.N+(1:G.ns)), ...
            'thetaW',a.thetaW,'thetaI',a.thetaI,'U',X(G.uIndex)); %#ok<AGROW>
        if p.verbose
            fprintf('day %8.2f: T=[%7.3f,%7.3f], h=[%8.2f,%8.2f], iceMax=%.4f, leak=%.4g m3/day/m, residual=%.2g\n', ...
                t/p.day,min(X(1:G.N)),max(X(1:G.N)),min(X(G.N+(1:G.ns))),max(X(G.N+(1:G.ns))),max(a.thetaI),D.leak*p.day,info.norm);
        end
        if p.saveOutput
            checkpoint=makeCheckpoint(X,t,G,p,initInfo,nrejected); %#ok<NASGU>
            checkpoint.runInventory=runInventory();
            writeCheckpointAtomic(p.outputPrefix,checkpoint,hist,cycles,snaps, ...
                leakHistSegment,leakHistJoint,leakHistBody,jointFreezeHist,jointHydraulicHist, ...
                mechanicalHistory,jointMechanicalHistory,strainHistory);
        end
        nextReport=nextReport+p.outputEveryDays*p.day;
    end
    if t>=nextYear-1e-7
        a=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
        ar=state(yearRef(1:G.N),yearRef(G.N+(1:G.ns)),G,p,yearRef(G.uIndex));
        eT=max(abs(X(1:G.N)-yearRef(1:G.N)));
        eH=max(abs(X(G.N+(1:G.ns))-yearRef(G.N+(1:G.ns))));
        eW=max(abs(a.total-ar.total)); eI=max(abs(a.thetaI-ar.thetaI));
        cycles=[cycles; (t-t0)/(p.day*p.yearDays),eT,eH,eW,eI]; %#ok<AGROW>
        meets=eT<p.spinup.tolT && eH<p.spinup.tolH && eW<p.spinup.tolWater && eI<p.spinup.tolIce;
        if size(cycles,1)>=p.spinup.consecutiveYears
            c=cycles(end-p.spinup.consecutiveYears+1:end,2:5);
            meets=meets && all(all(bsxfun(@lt,c,[p.spinup.tolT p.spinup.tolH p.spinup.tolWater p.spinup.tolIce])));
        else
            meets=false;
        end
        if p.verbose, fprintf('  annual same-phase changes: T %.4g C; h %.4g m; water %.4g; ice %.4g\n',eT,eH,eW,eI); end
        if strcmp(p.initial.mode,'spinup') && meets && ~cycleConverged
            cycleConverged=true;
            productionEnd=t+p.spinup.productionYears*p.yearDays*p.day;
            endTime=productionEnd;
            if p.verbose, fprintf('  periodic criterion met; now storing %g production years.\n',p.spinup.productionYears); end
        end
        yearRef=X; nextYear=nextYear+p.yearDays*p.day;
    end
end
a=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
out=struct('version','THM_v16_20260930_analytic_p1','T',X(1:G.N),'h',X(G.N+(1:G.ns)), ...
    'Psoil',p.rho*p.g*X(G.N+(1:G.ns)),'soilNodes',G.ids,'state',a,'pAll',G.xy,'mesh',S, ...
    'timeDay',t/p.day,'history',hist,'cycles',cycles,'snapshots',snaps,'parameters',p, ...
    'initialInfo',initInfo,'initialX',start,'periodicConverged',cycleConverged);
out.finalX=X;
out.algorithmRevision='THM_v16_research_20261002_local_cycle_coarsening';
out.soilConductivityField=G.soilConductivityField;
out.jointDamage=G.jointDamage;
out.jointCycleDamage=struct('enabled',G.jointCycleSettings.enabled, ...
    'settings',G.jointCycleSettings,'initialMemory',initInfo.jointCycleInitialMemory, ...
    'finalMemory',G.jointCycleMemory,'appliedSeverity',G.jointCycleAppliedSeverity, ...
    'history',jointCycleHistory,'sampleIceHistory',jointCycleIceHistory,'resolvedLaw',G.jointLaw, ...
    'historyBlocks',{{'day','appliedSeverity','queuedSeverity','meanCycleCount', ...
    'meanDose','connectedVolume_m2','referenceConductance_m_s'}}, ...
    'scope','Optional illustrative constant-volume channel coarsening; accepted-cycle lag, no calibrated gasket fatigue law.');
out.historyColumns={'day','dt_day','iterations','water_m2','energy_J_per_m','Qout_m2s', ...
    'Eout_W_per_m','waterError_m2','energyError_J_per_m','leak_m2s','scaledResidual', ...
    'jacobianBuilds','linearSolves','residualEvaluations'};
out.solverWork=struct('acceptedSteps',size(hist,1),'nonlinearIterations',sum(hist(:,3)), ...
    'jacobianBuilds',sum(hist(:,12)),'linearSolves',sum(hist(:,13)), ...
    'residualEvaluations',sum(hist(:,14)),'rejectedSteps',nrejected);
out.cycleColumns={'year','maxDeltaT_C','maxDeltaH_m','maxDeltaWater','maxDeltaIce'};
out.hydro=darcy(X,G,p);
out.numerics=struct('negativeHydraulicEdges',G.negativeHydraulicEdges, ...
    'negativeLiningEdges',G.negativeLiningEdges, ...
    'minimumGeometricEdgeWeight',G.minimumGeometricEdgeWeight, ...
    'advectionFormulation','Physical median-dual Darcy subface upwind', ...
    'phaseQuadratureOrder',p.integration.soilQuadratureOrder, ...
    'phaseQuadraturePoints',size(G.scv.P,1), ...
    'concreteHeatCapacity','Consistent Q4/T3 finite element mass matrix');
out.numerics.maximumAcceptedLocalWaterDefect_m2=maxAcceptedLocalWaterDefect;
out.numerics.maximumAcceptedLocalEnergyDefect_Jm=maxAcceptedLocalEnergyDefect;
out.leakage=struct('segmentRate_m2_per_s',leakHistSegment, ...
    'jointRate_m2_per_s',leakHistJoint,'bodyRate_m2_per_s',leakHistBody, ...
    'jointAnglesRad',p.leakage.jointAnglesRad(:)', ...
    'segmentCount',p.leakage.segmentCount,'interiorHead_m',p.leakage.interiorHead, ...
    'segmentBoundaryAnglesRad',B.segmentBoundaryAnglesRad(:)', ...
    'segmentPartitionMode',B.segmentPartitionMode);
out.leakage=summarize_segment_leakage(out.history,out.leakage,t0/p.day);
out.leakage.jointFreeze=freezeHistory(jointFreezeHist,p,G);
if G.jointCycleSettings.enabled
    out.leakage.jointFreeze.conductance_m_per_s= ...
        jointFreezeHist(:,1:8).*jointCycleHistory(:,42:49);
end
nj=numel(p.leakage.jointAnglesRad);
out.leakage.hydraulics=struct('history',jointHydraulicHist, ...
    'drivingHead_m',jointHydraulicHist(:,1:nj), ...
    'referenceHead_m',jointHydraulicHist(:,nj+(1:nj)), ...
    'unrestrictedRate_m2_per_s',jointHydraulicHist(:,2*nj+(1:nj)), ...
    'waterTableDepth_m',jointHydraulicHist(:,3*nj+1), ...
    'minTrans_m_per_s',jointHydraulicHist(:,3*nj+2), ...
    'negativeTransEdges',jointHydraulicHist(:,3*nj+3));
out.leakage.jointQuadrature=B.jointQuadrature;
if p.leakage.jointFreeze.enabled
    out.leakage.jointThermalPaths=B.jointThermalPaths;
    out.leakage.channelGeometry=G.jointChannelGeometry;
end
if p.thm.enabled
    out.mechanics=thmMechanicalOutput(X,G,p);
    out.mechanics.history=mechanicalHistory;
    out.mechanics.jointHistory=jointMechanicalHistory;
    out.mechanics.strainHistory=strainHistory;
    out.mechanics.strainHistoryColumns={'day','maxIncrementalVolumetricStrain', ...
        'maxIncrementalPrincipalStrain','maxIncrementalTensorShearStrain', ...
        'maxIncrementalGradientNorm','maxAbsolutePrincipalStrain','maxAbsoluteGradientNorm'};
    out.mechanics.jointHistoryBlocks={'day','minGap_m','meanGap_m','maxGap_m', ...
        'maxAbsSlip_m','maxContactPressure_Pa','boltForce_Nm', ...
        'capillaryConductance_m_s','openSlitConductance_m_s','totalConductance_m_s'};
    out.mechanics.historyColumns={'day','maxAbsVolumetricStrain','maxJointGap_m','maxContactPressure_Pa','maxFreeMechanicalResidual_Nm','maxOmittedSlitVolumeFraction'};
elseif p.mechanics.enable
    out.mechanics=mechanics(X,start,G,S,p);
end
if strcmp(p.initial.mode,'spinup') && ~cycleConverged
    warning('TH:spinup','Maximum spin-up duration reached WITHOUT meeting periodic tolerances. Do not label the state independent of initialization.');
end
if p.saveOutput
    save([p.outputPrefix '_results.mat'],'out','-v7');
    checkpoint=makeCheckpoint(X,t,G,p,initInfo,nrejected); %#ok<NASGU>
    checkpoint.runInventory=runInventory();
    writeCheckpointAtomic(p.outputPrefix,checkpoint,hist,cycles,snaps, ...
        leakHistSegment,leakHistJoint,leakHistBody,jointFreezeHist,jointHydraulicHist, ...
        mechanicalHistory,jointMechanicalHistory,strainHistory);
end
if p.plot, plotOutput(out,G); end
    function run=runInventory()
        reportTime=nextReport;
        while reportTime<=t+1e-7, reportTime=reportTime+p.outputEveryDays*p.day; end
        yearTime=nextYear; ref=yearRef;
        if yearTime<=t+1e-7, yearTime=yearTime+p.yearDays*p.day; ref=X; end
        run=struct('startX',start,'startTime',t0,'endTime',endTime, ...
            'nextReport',reportTime,'nextYear',yearTime,'yearRef',ref, ...
            'nextDt',dt,'priorStart',priorStart,'priorAcceptedDt',priorAcceptedDt, ...
            'maxLocalWaterDefect',maxAcceptedLocalWaterDefect, ...
            'maxLocalEnergyDefect',maxAcceptedLocalEnergyDefect);
        if G.jointCycleSettings.enabled
            run.jointCycleHistory=jointCycleHistory;
            run.jointCycleIceHistory=jointCycleIceHistory;
            run.startCycleMemory=initInfo.jointCycleInitialMemory;
        end
    end
end

function nextDt=acceptedTimeStep(proposedDt,acceptedDt,info,p)
% Count Newton iterations, not mesh-dependent finite-difference evaluations.
% An accepted step below the high-iteration threshold always regrows; a
% truncated monthly/yearly remainder cannot reduce the requested step.
iterationFraction=max(info.iter,1)/p.solver.maxIter;
if iterationFraction>=p.solver.dtHighIterationFraction
    factor=.7;
else
    factor=p.solver.dtGrowthSafety*sqrt( ...
        p.solver.dtTargetNewtonFraction/iterationFraction);
    factor=min(p.solver.dtGrowthMax,max(p.solver.dtAcceptedGrowthFloor,factor));
end
if acceptedDt<proposedDt*(1-1e-10), factor=max(1,factor); end
nextDt=min(p.dtMax,max(p.dtMin,proposedDt*factor));
end

function writeCheckpointAtomic(prefix,checkpoint,hist,cycles,snaps, ...
    leakHistSegment,leakHistJoint,leakHistBody,jointFreezeHist,jointHydraulicHist, ...
    mechanicalHistory,jointMechanicalHistory,strainHistory)
% Keep the previous complete month available until the new MAT is written.
temporaryFile=[prefix '_checkpoint_tmp.mat'];
checkpointFile=[prefix '_checkpoint.mat'];
save(temporaryFile,'checkpoint','hist','cycles','snaps', ...
    'leakHistSegment','leakHistJoint','leakHistBody','jointFreezeHist','jointHydraulicHist', ...
    'mechanicalHistory','jointMechanicalHistory','strainHistory','-v7');
[ok,message]=movefile(temporaryFile,checkpointFile,'f');
assert(ok,message);
end

function p=defaults()
p.mesh=[]; p.meshFile=fullfile(fileparts(mfilename('fullpath')),'example_mesh_repaired_interface.mat');
p.day=86400; p.yearDays=365; p.durationDays=730;
p.dtInitial=.1*p.day; p.dtMax=p.day; p.dtMin=p.day/4096;
p.outputEveryDays=30; p.outputPrefix='TH_v16_cvfe_phase'; p.saveOutput=true; p.plot=false; p.verbose=true;
p.rho=1000; p.g=9.81; p.L=3.34e5; p.Tm=273.15;
p.thm=thmDefaults();
p.enthalpyReference_Jkg=0; % arbitrary common water/ice energy datum, for invariance verification
p.integration=struct('soilQuadratureOrder',5,'subdivisionLevel',0);
p.porosity=.35; p.thetaR=.05; p.vgAlpha=.65; p.vgN=1.67; % alpha in 1/m, illustrative silt
p.Ksat=1e-14/1e-3*p.rho*p.g; % m/s, inherited intrinsic permeability
p.soilKsatField=struct('coordinates_m',[],'logMultiplier',[], ...
    'quadratureCoordinates_m',[],'quadratureLogMultiplier',[]);
p.soilViscosity=struct('enabled',true,'referenceTemperature_C',20);
p.Ss=1e-5; % 1/m, saturated elastic storage; SAME storage in mass and enthalpy
p.iceImpedance=0; % extra empirical ice blockage only if measured and calibrated
p.cw=4180; p.ci=2100; p.rhoS=2650; p.cs=800;
p.kw=.57; p.ki=2.2; p.ks=2.3; p.ka=.025;
p.rhoC=2500; p.cc=900; p.kc=1.7;
p.topMin=-10; p.topMax=15; p.tunnelMin=-5; p.tunnelMax=8;
p.topHTC=15; p.tunnelHTC=10; % air-to-surface Robin W/m2/K (Inf -> Dirichlet)
p.geothermalFlux=2.3*.006; % independently prescribed basal heat input W/m2
p.waterTableDepth=5; p.hydroFar=true; p.topSeepage=true;
p.forcing=struct('day',[],'topAirC',[],'tunnelAirC',[],'waterTableDepth_m',[]);
p.groundwaterInflowT=2.5; p.groundwaterInflowGradient=.006; % independent of initial T
p.topLeakConductance=1e-6; % m/s per m pressure head; one-way regularized seepage
p.tunnelDrain=true; p.liningThickness=.45;
p.liningK=1e-10; % radial lining hydraulic conductivity m/s; illustrative
p.leakage=struct('segmentCount',6,'segmentStartAngleRad',0, ...
    'interiorHead',0, ...
    'jointAnglesRad',[],'jointWidth_m',[],'jointConductance',[],'jointEntryHead',[], ...
    'jointFreeze',joint_freeze_defaults(), ...
    'jointStorage',struct('enabled',true,'tortuosity',1), ...
    'jointDamage',struct('enabled',false,'severity',0, ...
        'radiusGrowthAtUnitDamage',.5,'channelDensityGrowthAtUnitDamage',0), ...
    'jointCycleDamage',joint_cycle_damage_defaults());
p.interfaceHTC=5;
p.advection=true;
p.initial=struct('mode','profile','surfaceT',-1,'gradient',.006,'T',[],'h',[], ...
    'checkpointFile','','startDay',0,'resumeHistory',false, ...
    'allowNumericalReintegration',false);
p.spinup=struct('tolT',.02,'tolH',.02,'tolWater',1e-4,'tolIce',1e-4, ...
    'consecutiveYears',2,'productionYears',1);
p.solver=struct('stepU',.001,'maxIter',40,'tol',1e-8,'stepT',2,'stepH',50,'fd',1e-8, ...
    'lineSearch',16,'massScale',.01,'energyScale',1e7, ...
    'jacobianReuse',4,'jacobianRefreshRatio',.65,'predictor',true, ...
    'reuseAcrossSteps',false,'jacobian','analytic','reuseAcceptedResidual',true, ...
    'reuseConstitutiveState',true, ... % exact old/accepted SCV and contact state cache
    'dtMaxByMonth_days',[], ... % optional calendar caps; no effect when empty
    'scvTangentAssembly','element', ... % same quadrature tangent, fewer sparse products
    'dtTargetResidualEvals',180, ... % legacy input; no longer controls the accepted timestep
    'dtTargetNewtonFraction',.4,'dtAcceptedGrowthFloor',1.05, ...
    'dtHighIterationFraction',.75,'dtGrowthSafety',1.10, ...
    'dtGrowthMax',1.8,'dtShrinkMin',.5); % dtShrinkMin is also legacy
p.mechanics=struct('enable',true,'Esoil',80e6,'nuSoil',.35,'Econc',34.5e9, ...
    'nuConc',.2,'alphaSoil',1.2e-5,'alphaConc',1e-5,'biot',.65, ...
    'kn',5e9,'kt',1e9,'frostEfficiency',0);
end

function a=merge(a,b)
f=fieldnames(b);
for i=1:numel(f)
    k=f{i};
    if ~isfield(a,k), error('Unknown option: %s',k); end
    if isstruct(a.(k)) && isstruct(b.(k)), a.(k)=merge(a.(k),b.(k)); else, a.(k)=b.(k); end
end
end

function validateParameters(p)
assert(isscalar(p.initial.allowNumericalReintegration) ...
    &&ismember(p.initial.allowNumericalReintegration,[false true]), ...
    'THM:ReintegrationOption','allowNumericalReintegration must be logical 0/1.');
assert(~p.initial.allowNumericalReintegration|| ...
    (strcmp(p.initial.mode,'checkpoint')&&~p.initial.resumeHistory), ...
    'THM:ReintegrationHistory', ...
    'Reintegration requires a new checkpoint-seeded run; saved history cannot be reused.');
assert(p.thetaR>=0 && p.thetaR<p.porosity && p.porosity<1);
assert(p.vgN>1 && p.vgAlpha>0 && p.Ksat>0 && p.Ss>0);
assert(p.dtMin>0 && p.dtInitial>=p.dtMin && p.dtMax>=p.dtInitial);
assert(p.outputEveryDays>0 && p.durationDays>0);
assert(p.topHTC>0 && p.tunnelHTC>0 && p.interfaceHTC>0);
assert(p.initial.gradient>=0 && p.rho>0 && p.L>0);
assert(ismember(p.integration.soilQuadratureOrder,[2 5]));
assert(isscalar(p.integration.subdivisionLevel)&&p.integration.subdivisionLevel>=0 ...
    &&p.integration.subdivisionLevel<=2&&p.integration.subdivisionLevel==fix(p.integration.subdivisionLevel));
assert(isscalar(p.enthalpyReference_Jkg)&&isfinite(p.enthalpyReference_Jkg));
if p.thm.enabled
    assert(p.leakage.jointFreeze.enabled&&p.leakage.jointStorage.enabled, ...
        'THM requires capillary thermal paths and inventory (iceEnabled may be false).');
    assert(p.mechanics.biot>=0&&p.mechanics.biot<=1);
    assert(p.mechanics.frostEfficiency==0,'Do not double-count a prescribed frost eigenstrain.');
    assert(p.thm.openingTransmissionFraction>=0&&p.thm.openingTransmissionFraction<=1);
    assert(p.thm.maximumVolumetricStrain>0&&p.thm.maximumVolumetricStrain<=.05);
    assert(p.thm.maximumPrincipalStrain>0&&p.thm.maximumPrincipalStrain<=.05);
    assert(p.thm.maximumDisplacementGradient>0&&p.thm.maximumDisplacementGradient<=.1);
    assert(p.thm.maxOmittedSlitVolumeFraction>=0&&p.thm.maxOmittedSlitVolumeFraction<=.01);
    assert(~strcmp(p.initial.mode,'spinup'),'THM:SpinupUnsupported', ...
        'TH-only periodic shortcut is not validated for contact THM; use run_TH_segment_year20_thm for full history.');
end
assert(any(strcmp(p.initial.mode,{'profile','fields','checkpoint','spinup'})));
assert(p.liningThickness>0 && p.liningK>=0 && p.iceImpedance>=0);
assert(p.solver.jacobianReuse>=0 && p.solver.jacobianReuse==fix(p.solver.jacobianReuse));
assert((islogical(p.solver.reuseAcceptedResidual)||isnumeric(p.solver.reuseAcceptedResidual)) ...
    &&isscalar(p.solver.reuseAcceptedResidual)&&isreal(p.solver.reuseAcceptedResidual) ...
    &&isfinite(p.solver.reuseAcceptedResidual)&&ismember(p.solver.reuseAcceptedResidual,[0 1]), ...
    'THM:ResidualReuse','reuseAcceptedResidual must be scalar logical 0/1.');
assert(p.solver.jacobianRefreshRatio>0 && p.solver.jacobianRefreshRatio<1);
assert(isscalar(p.solver.reuseConstitutiveState)&& ...
    ismember(p.solver.reuseConstitutiveState,[0 1]), ...
    'THM:StateReuse','reuseConstitutiveState must be scalar logical 0/1.');
caps=p.solver.dtMaxByMonth_days;
assert(isnumeric(caps)&&isreal(caps)&&(isempty(caps) ...
    ||(isvector(caps)&&numel(caps)==12&&all(isfinite(caps)) ...
    &&all(caps*p.day>=p.dtMin)&&all(caps*p.day<=p.dtMax))), ...
    'THM:CalendarStepCap','Supply empty caps or 12 positive monthly day caps within dtMin/dtMax.');
assert(any(strcmp(p.solver.jacobian,{'analytic','coloredFD'})), ...
    'THM:JacobianMode','Use analytic or coloredFD Jacobian.');
assert(any(strcmp(p.solver.scvTangentAssembly,{'element','sparse'})), ...
    'THM:SCVTangentAssembly','Use element or sparse SCV tangent assembly.');
assert(p.solver.dtTargetResidualEvals>0 && p.solver.dtGrowthSafety>0);
assert(p.solver.dtTargetNewtonFraction>0 && p.solver.dtTargetNewtonFraction<1);
assert(p.solver.dtAcceptedGrowthFloor>=1 && ...
    p.solver.dtAcceptedGrowthFloor<=p.solver.dtGrowthMax);
assert(p.solver.dtHighIterationFraction>p.solver.dtTargetNewtonFraction && ...
    p.solver.dtHighIterationFraction<=1);
assert(p.solver.dtGrowthMax>=1 && p.solver.dtShrinkMin>0 && p.solver.dtShrinkMin<=1);
L=p.leakage; nj=numel(L.jointAnglesRad);
assert(isscalar(L.segmentCount)&&L.segmentCount>=1&&L.segmentCount==fix(L.segmentCount));
assert(isscalar(L.segmentStartAngleRad)&&isfinite(L.segmentStartAngleRad));
assert(isscalar(L.interiorHead)&&isfinite(L.interiorHead));
assert(all(isfinite(L.jointAnglesRad(:)))&&isreal(L.jointAnglesRad));
assert(all(isfinite(L.jointConductance(:)))&&all(L.jointConductance(:)>=0));
assert(all(isfinite(L.jointWidth_m(:)))&&all(L.jointWidth_m(:)>0));
assert(all(isfinite(L.jointEntryHead(:)))&&all(L.jointEntryHead(:)>=0));
assert(nj==0 || numel(L.jointWidth_m)==1 || numel(L.jointWidth_m)==nj, ...
    'Supply one positive physical width per joint or a shared scalar.');
assert(nj==0 || (numel(L.jointConductance)==1 || numel(L.jointConductance)==nj), ...
    'Supply one conductance per joint or a shared scalar.');
assert(nj==0 || isempty(L.jointEntryHead) || numel(L.jointEntryHead)==1 || ...
    numel(L.jointEntryHead)==nj,'Joint entry heads must be scalar or one per joint.');
assert(nj>0 || isempty(L.jointConductance),'Joint conductance requires joint angles.');
assert(nj>0 || isempty(L.jointWidth_m),'Joint width requires joint angles.');
if L.jointFreeze.enabled
    assert(nj>0,'Joint freezing requires at least one joint.');
    joint_freeze_conductance(zeros(nj,L.jointFreeze.quadraturePoints),L.jointFreeze);
end
assert(islogical(p.soilViscosity.enabled)||ismember(p.soilViscosity.enabled,[0 1]));
water_viscosity_iapws(p.soilViscosity.referenceTemperature_C);
assert(islogical(L.jointStorage.enabled)||ismember(L.jointStorage.enabled,[0 1]));
assert(all(isfinite(L.jointStorage.tortuosity(:)))&&all(L.jointStorage.tortuosity(:)>=1));
forcing=p.forcing;
series={'topAirC','tunnelAirC','waterTableDepth_m'};
hasSeries=false;
for k=1:numel(series)
    values=forcing.(series{k});
    if isempty(values), continue; end
    hasSeries=true;
    assert(isvector(values)&&numel(values)==numel(forcing.day)&&all(isfinite(values(:))), ...
        'Each observed forcing must match forcing.day.');
    assert(abs(values(1)-values(end))<1e-8, ...
        'Repeated annual forcing must have identical first and last values.');
end
if hasSeries
    assert(isvector(forcing.day)&&numel(forcing.day)>=2&& ...
        abs(forcing.day(1))<1e-10&&abs(forcing.day(end)-p.yearDays)<1e-10&& ...
        all(diff(forcing.day)>0),'forcing.day must increase from 0 to yearDays.');
end
end

function G=geometry(S,p)
[G.jointLaw,G.jointDamage]=joint_microchannel_damage(p.leakage);
G.fixedJointLaw=G.jointLaw;
G.jointCycleSettings=joint_cycle_damage_defaults(p.leakage.jointCycleDamage);
G.jointCycleMemory=[];G.jointCycleAppliedSeverity=[];
if G.jointCycleSettings.enabled
    assert(p.thm.enabled&&p.advection&&p.leakage.jointFreeze.enabled ...
        &&p.leakage.jointFreeze.iceEnabled&&p.leakage.jointStorage.enabled, ...
        'JointCycle:Physics','Coarsening requires THM, advection, freezing and channel inventory.');
    assert(~G.jointDamage.enabled,'JointCycle:Physics', ...
        'Use initialSeverity for evolving cases; fixed damage cannot be stacked ambiguously.');
    assert(~strcmp(p.initial.mode,'spinup'),'JointCycle:Physics', ...
        'Irreversible coarsening cannot use periodic spinup termination.');
    G=applyJointCycleGeometry(G,G.jointCycleSettings.initialSeverity);
end
G.xy=S.pAll; G.N=size(G.xy,1);
if isfield(S,'xc')&&isfield(S,'yc')
    G.center=[S.xc S.yc];
else
    assert(~isfield(S,'soilTunnelCoupleTable')||isempty(S.soilTunnelCoupleTable), ...
        'Segmented lining mesh must supply tunnel center xc and yc.');
    G.center=(min(G.xy,[],1)+max(G.xy,[],1))/2;
end
if isfield(S,'coordInfo') && isfield(S.coordInfo,'sectionColumns')
    if isempty(strfind(lower(S.coordInfo.sectionColumns),'down')) %#ok<STREMP>
        error('Mesh must explicitly use Z positive down. Regenerate with supplied unified mesh script.');
    end
end
G.tri=S.tSoil; G.ids=unique(G.tri(:)); G.ns=numel(G.ids);
G.map=zeros(G.N,1); G.map(G.ids)=1:G.ns; G.local=G.map(G.tri);
[G.area,G.grad]=triGeometry(G.xy,G.tri);
G.vol=accumarray(G.local(:),repmat(G.area/3,3,1),[G.ns 1]);
G.volAll=zeros(G.N,1); G.volAll(G.ids)=G.vol;
[G.ri,G.ci,G.kvals]=triTriplets(G.tri,G.area,G.grad);
G.Kconc=sparse(G.N,G.N); G.Mconc=sparse(G.N,G.N);
G.vconc=zeros(G.N,1); G.qsets={}; G.ctri=zeros(0,3);
for name={'qTunnelGlobal','qInvertGlobal'}
    if isfield(S,name{1}) && ~isempty(S.(name{1}))
        qs=S.(name{1}); G.qsets{end+1}=qs;
        for e=1:size(qs,1)
            ids=qs(e,:); [K,v,M]=quadThermal(G.xy(ids,:));
            G.Kconc(ids,ids)=G.Kconc(ids,ids)+p.kc*K;
            G.Mconc(ids,ids)=G.Mconc(ids,ids)+M;
            G.vconc(ids)=G.vconc(ids)+v;
        end
    end
end
if isfield(S,'tInvertTriGlobal') && ~isempty(S.tInvertTriGlobal)
    G.ctri=S.tInvertTriGlobal; [ar,gr]=triGeometry(G.xy,G.ctri);
    [ii,jj,vv]=triTriplets(G.ctri,ar,gr);
    G.Kconc=G.Kconc+sparse(ii,jj,p.kc*vv,G.N,G.N);
    for e=1:size(G.ctri,1)
        ids=G.ctri(e,:);
        G.Mconc(ids,ids)=G.Mconc(ids,ids)+ar(e)/12*[2 1 1;1 2 1;1 1 2];
    end
    G.vconc=G.vconc+accumarray(G.ctri(:),repmat(ar/3,3,1),[G.N 1]);
end
assert(all(G.volAll+G.vconc>0),'Unused nodes or zero-volume elements in mesh.');
G.C=G.Kconc; G.interface=zeros(0,4);
if isfield(S,'soilTunnelCoupleTable') && ~isempty(S.soilTunnelCoupleTable)
    C=S.soilTunnelCoupleTable; G.interface=C(:,[2 3 6 7]);
    for e=1:size(G.interface,1)
        ids=G.interface(e,:); xs=G.xy(ids(1:2),:); xt=G.xy(ids(3:4),:);
        if norm(xs-xt,'fro')>1e-8
            if norm(xs-flipud(xt),'fro')<1e-8
                ids(3:4)=ids([4 3]); G.interface(e,:)=ids;
            else
                error('Soil-lining edge mismatch at interface row %d.',e);
            end
        end
        L=norm(xs(2,:)-xs(1,:)); M=L/6*[2 1;1 2];
        G.C(ids,ids)=G.C(ids,ids)+p.interfaceHTC*[M -M;-M M];
    end
end
% Global pairwise FEM hydraulic flux. Each pair carries one antisymmetric flux.
pairs=[G.local(:,[1 2]);G.local(:,[2 3]);G.local(:,[3 1])];
vals=-[G.area.*sum(G.grad(:,:,1).*G.grad(:,:,2),2); ...
       G.area.*sum(G.grad(:,:,2).*G.grad(:,:,3),2); ...
       G.area.*sum(G.grad(:,:,3).*G.grad(:,:,1),2)];
[G.edges,~,G.edgeMap]=unique(sort(pairs,2),'rows');
G.edgeBase=vals; G.edgeElement=repmat((1:numel(G.area))',3,1);
geometricWeight=accumarray(G.edgeMap,G.edgeBase,[size(G.edges,1) 1]);
G.negativeHydraulicEdges=sum(geometricWeight< -1e-12);
G.minimumGeometricEdgeWeight=min(geometricWeight);
liningPairs=sort(G.map(G.interface(:,1:2)),2);
if isempty(liningPairs)
    G.negativeLiningEdges=0;
else
    G.negativeLiningEdges=sum(geometricWeight< -1e-12 & ...
        ismember(G.edges,liningPairs,'rows'));
end
G.inc=sparse([G.edges(:,1);G.edges(:,2)],[(1:size(G.edges,1))';(1:size(G.edges,1))'], ...
    [ones(size(G.edges,1),1);-ones(size(G.edges,1),1)],G.ns,size(G.edges,1));
% Exact linear map from nodal conductivity to each edge transmissivity.
% It replaces repeated element averaging and ACCUMARRAY during Newton.
ne=size(G.tri,1); nf=size(G.edges,1);
triangleMean=sparse(repmat((1:ne)',3,1),G.local(:),1/3,ne,G.ns);
edgeToTriangle=sparse(G.edgeMap,G.edgeElement,G.edgeBase,nf,ne);
G.edgeKPattern=edgeToTriangle*triangleMean;
G.edgeElementPattern=edgeToTriangle;
G.scv=build_soil_scv_quadrature(G.xy(G.ids,:),G.local, ...
    p.integration.soilQuadratureOrder,p.integration.subdivisionLevel);
G.scv.tangentAssembly=build_scv_tangent_assembly(G.scv,G.local);
G.cvfe=build_soil_cvfe_flux(G.xy,G.tri,G.ids);
G.soilSelect=sparse((1:G.ns)',G.ids,1,G.ns,G.N);
G.soilConductivityField=resolve_soil_conductivity_field( ...
    p.soilKsatField,G.xy(G.ids,:),G.scv.P);
A=sparse(G.ri,G.ci,1,G.N,G.N)+spones(G.C)+speye(G.N);
A=spones(A+A'); H=A(G.ids,G.ids);
% Boundary-water enthalpy reactions only couple immediate neighboring soil nodes.
pattern=[A,A(:,G.ids);A(G.ids,:),H]; G.pattern=spones(pattern);
G.groups=colorColumns(G.pattern);
[G.jrow,G.jcol]=find(G.pattern);
% Precompute compressed-FD scatter maps once.  The v5 implementation called
% ISMEMBER for every colour and every Newton build, which is pure overhead.
G.jnz=cell(size(G.groups)); G.jloc=cell(size(G.groups));
for c=1:numel(G.groups)
    [tf,loc]=ismember(G.jcol,G.groups{c});
    G.jnz{c}=find(tf); G.jloc{c}=loc(tf);
end
G.rs=[(G.volAll+G.vconc)*p.solver.energyScale;G.vol*p.solver.massScale];
G.xs=[ones(G.N,1);10*ones(G.ns,1)];
G.z0=min(G.xy(:,2)); G.zWT=G.z0+p.waterTableDepth;
G.uIndex=zeros(0,1); G.nm=0;
end

function G=jointCouplingPattern(G,B)
% Remote lining temperatures affect joint mass and drain-enthalpy rows.
P=B.jointThermalPaths.interpolation;
nj=B.jointThermalPaths.nj; nq=B.jointThermalPaths.nq;
jointNodes=spones(kron(speye(nj),ones(1,nq))*spones(P));
QT=spones(spones(B.jointNodeWeight')*jointNodes);
extra=[G.soilSelect'*QT,sparse(G.N,G.ns);QT,sparse(G.ns,G.ns)];
thermal=spones(jointNodes'*jointNodes);
thermalSoil=spones(jointNodes'*spones(B.jointNodeWeight));
extra=extra+[thermal+thermalSoil*G.soilSelect,thermalSoil; ...
    sparse(G.ns,G.N+G.ns)];
G.pattern=spones(G.pattern+extra);
G.groups=colorColumns(G.pattern); [G.jrow,G.jcol]=find(G.pattern);
G.jnz=cell(size(G.groups)); G.jloc=cell(size(G.groups));
for c=1:numel(G.groups)
    [tf,loc]=ismember(G.jcol,G.groups{c});
    G.jnz{c}=find(tf); G.jloc{c}=loc(tf);
end
end

function f=freezeHistory(history,p,G)
nj=numel(p.leakage.jointAnglesRad);
c=G.jointLaw.jointConductance(:)';
if nj==0, c=zeros(1,0); elseif isscalar(c), c=repmat(c,1,nj); end
f=struct('enabled',p.leakage.jointFreeze.enabled,'history',history, ...
    'factor',history(:,1:nj),'minTemperature_C',history(:,nj+(1:nj)), ...
    'maxTemperature_C',history(:,2*nj+(1:nj)), ...
    'frozenLengthFraction',history(:,3*nj+(1:nj)), ...
    'conductance_m_per_s',history(:,1:nj).*c);
end

function [ar,gr]=triGeometry(x,t)
a=x(t(:,1),:); b=x(t(:,2),:); c=x(t(:,3),:);
a2=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(c(:,1)-a(:,1)).*(b(:,2)-a(:,2));
assert(all(abs(a2)>1e-12),'Degenerate soil/structure triangle.');
ar=abs(a2)/2; gr=zeros(size(t,1),2,3);
% N1 gradient = [z2-z3, y3-y2]/A2
gr(:,:,1)=[b(:,2)-c(:,2), c(:,1)-b(:,1)]./a2;
gr(:,:,2)=[c(:,2)-a(:,2), a(:,1)-c(:,1)]./a2;
gr(:,:,3)=[a(:,2)-b(:,2), b(:,1)-a(:,1)]./a2;
end

function [ii,jj,v]=triTriplets(t,a,g)
ne=size(t,1); ii=zeros(ne,9); jj=ii; v=ii; k=0;
for i=1:3
    for j=1:3
        k=k+1; ii(:,k)=t(:,i); jj(:,k)=t(:,j); v(:,k)=a.*sum(g(:,:,i).*g(:,:,j),2);
    end
end
ii=ii(:); jj=jj(:); v=v(:);
end

function [K,v,M]=quadThermal(x)
K=zeros(4); v=zeros(4,1); M=zeros(4); signs=[];
for xi=[-1 1]/sqrt(3)
    for eta=[-1 1]/sqrt(3)
        [N,d]=qshape(xi,eta); J=d*x; dj=det(J); signs(end+1)=sign(dj); %#ok<AGROW>
        assert(abs(dj)>1e-12,'Degenerate Q4.'); dn=(J\d)';
        K=K+dn*dn'*abs(dj); v=v+N'*abs(dj); M=M+N'*N*abs(dj);
    end
end
assert(all(signs==signs(1)),'Folded Q4: Jacobian changes sign.');
end

function [N,d]=qshape(xi,eta)
N=.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
d=.25*[-(1-eta),1-eta,1+eta,-(1+eta);-(1-xi),-(1+xi),1+xi,1-xi];
end

function groups=colorColumns(P)
C=spones(P'*P); n=size(P,2); col=zeros(n,1);
[~,order]=sort(full(sum(C,1)),'descend');
for j=order
    used=col(find(C(:,j))); used=used(used>0);
    k=1; while any(used==k), k=k+1; end
    col(j)=k;
end
groups=cell(max(col),1);
for k=1:numel(groups), groups{k}=find(col==k); end
end

function B=boundaries(G,S,p)
x=G.xy; tol=1e-8; ids=G.ids;
B.top=ids(abs(x(ids,2)-min(x(ids,2)))<tol);
B.bottom=ids(abs(x(ids,2)-max(x(ids,2)))<tol);
B.left=ids(abs(x(ids,1)-min(x(ids,1)))<tol);
B.right=ids(abs(x(ids,1)-max(x(ids,1)))<tol);
B.topW=lineWeights(x,B.top,1); B.bottomW=lineWeights(x,B.bottom,1);
B.tunnelW=zeros(G.N,1); B.drainW=zeros(G.ns,1);
if isfield(S,'tunnelInnerRingGlobal') && ~isempty(S.tunnelInnerRingGlobal)
    ring=S.tunnelInnerRingGlobal(:); exposed=true(numel(ring),1);
    if isfield(S,'invertBottomRingGlobal') && ~isempty(S.invertBottomRingGlobal)
        covered=sort([S.invertBottomRingGlobal(1:end-1),S.invertBottomRingGlobal(2:end)],2);
        edges=[ring,ring([2:end 1])]; exposed=~ismember(sort(edges,2),covered,'rows');
    else
        edges=[ring,ring([2:end 1])];
    end
    B.tunnelW=edgeWeights(x,edges(exposed,:));
    if isfield(S,'bInvertTopGlobal') && ~isempty(S.bInvertTopGlobal)
        B.tunnelW=B.tunnelW+edgeWeights(x,S.bInvertTopGlobal);
    end
end
for e=1:size(G.interface,1)
    k=G.interface(e,1:2); L=norm(x(k(2),:)-x(k(1),:));
    B.drainW(G.map(k))=B.drainW(G.map(k))+L/2;
end
% Hydraulic head is H=h-z because the mesh uses positive-down elevation.
% The soil-to-tunnel path descends at the crown and rises at the invert.
% Check the specified lining thickness against the actual circular mesh.
B.liningZDrop=zeros(G.ns,1);
B.bodyRadialFactor=zeros(G.ns,1);
if ~isempty(G.interface)
    interfaceNodes=unique(G.interface(:,1:2));
    outerRadius=hypot(x(interfaceNodes,1)-G.center(1), ...
        x(interfaceNodes,2)-G.center(2));
    if isfield(S,'tunnelInnerRingGlobal') && ~isempty(S.tunnelInnerRingGlobal)
        innerNodes=S.tunnelInnerRingGlobal(:);
        innerRadius=hypot(x(innerNodes,1)-G.center(1), ...
            x(innerNodes,2)-G.center(2));
        assert(max(innerRadius)-min(innerRadius)<1e-6, ...
            'Total-head lining leakage requires a circular inner ring.');
        meshThickness=outerRadius-mean(innerRadius);
        assert(all(abs(meshThickness-p.liningThickness)<max(1e-4,.01*p.liningThickness)), ...
            'liningThickness must match the actual soil-to-tunnel radial thickness.');
    else
        assert(all(outerRadius>p.liningThickness), ...
            'Lining thickness exceeds the outer interface radius.');
    end
    B.liningZDrop(G.map(interfaceNodes))= ...
        -p.liningThickness*(x(interfaceNodes,2)-G.center(2))./outerRadius;
    % Exact steady radial conductance per unit outer arc length of an
    % annulus, with the input K interpreted as a material conductivity.
    innerPathRadius=outerRadius-p.liningThickness;
    assert(all(innerPathRadius>0),'Invalid annular lining thickness.');
    B.bodyRadialFactor(G.map(interfaceNodes))= ...
        1./(outerRadius.*log(outerRadius./innerPathRadius));
end
% Split each interface edge at segment boundaries. If one joint is supplied
% per segment, the actual joint angles define the segment boundaries;
% otherwise retain the legacy equal-angle statistical sectors.
B.sectorNodeWeight=sparse(G.ns,p.leakage.segmentCount);
if numel(p.leakage.jointAnglesRad)==p.leakage.segmentCount
    sectorBounds=sort(mod(p.leakage.jointAnglesRad(:)- ...
        p.leakage.segmentStartAngleRad,2*pi));
    assert(min(diff([sectorBounds;sectorBounds(1)+2*pi]))>1e-9, ...
        'Distinct segment joints are required to define physical segments.');
    B.segmentPartitionMode='joint_boundaries';
else
    sectorBounds=(0:p.leakage.segmentCount-1)'*2*pi/p.leakage.segmentCount;
    B.segmentPartitionMode='equal_sectors';
end
B.segmentBoundaryAnglesRad=mod(sectorBounds+p.leakage.segmentStartAngleRad,2*pi);
for edge=1:size(G.interface,1)
    nodes=G.interface(edge,1:2); local=G.map(nodes);
    theta=mod(atan2(x(nodes,2)-G.center(2),x(nodes,1)-G.center(1)) ...
        -p.leakage.segmentStartAngleRad,2*pi);
    if theta(2)-theta(1)>pi, theta(2)=theta(2)-2*pi; end
    if theta(2)-theta(1)<-pi, theta(2)=theta(2)+2*pi; end
    assert(abs(diff(theta))>1e-12,'Degenerate angular interface edge.');
    low=min(theta); high=max(theta);
    crossings=segmentCrossings(low,high,sectorBounds);
    fractions=sort([0;(crossings(:)-theta(1))/diff(theta);1]);
    len=norm(x(nodes(2),:)-x(nodes(1),:));
    for piece=1:numel(fractions)-1
        u=fractions(piece); v=fractions(piece+1); mid=(u+v)/2;
        sector=segmentIndex(theta(1)+(u+v)/2*diff(theta),sectorBounds);
        B.sectorNodeWeight(local,sector)=B.sectorNodeWeight(local,sector)+ ...
            len*(v-u)*[1-mid;mid];
    end
end
assert(max(abs(sum(B.sectorNodeWeight,2)-B.drainW))<1e-8, ...
    'Segment weights do not partition the lining interface.');
nj=numel(p.leakage.jointAnglesRad);
B.jointNodeWeight=sparse(nj,G.ns);
B.jointSectorNodeWeight=cell(p.leakage.segmentCount,1);
for sector=1:p.leakage.segmentCount
    B.jointSectorNodeWeight{sector}=sparse(nj,G.ns);
end
if nj>0
    assert(~isempty(G.interface),'Joint leakage requires a soil-lining interface.');
    widths=p.leakage.jointWidth_m(:);
    if isscalar(widths), widths=repmat(widths,nj,1); end
    edges=G.interface(:,1:2);
    radius=mean(hypot(x(edges(:,1),1)-G.center(1),x(edges(:,1),2)-G.center(2)));
    for firstJoint=1:nj
        for secondJoint=firstJoint+1:nj
            angularDistance=abs(mod(p.leakage.jointAnglesRad(firstJoint)- ...
                p.leakage.jointAnglesRad(secondJoint)+pi,2*pi)-pi);
            assert(radius*angularDistance>=(widths(firstJoint)+widths(secondJoint))/2-1e-10, ...
                'Physical joint bands overlap; combine them into one conductance.');
        end
    end
    for jnt=1:nj
        halfAngle=widths(jnt)/(2*radius);
        assert(halfAngle<pi/2,'Joint band cannot cover half the circumference.');
        centerAngle=mod(p.leakage.jointAnglesRad(jnt),2*pi);
        for edge=1:size(edges,1)
            nodes=edges(edge,:); local=G.map(nodes);
            theta=mod(atan2(x(nodes,2)-G.center(2),x(nodes,1)-G.center(1)),2*pi);
            if theta(2)-theta(1)>pi, theta(2)=theta(2)-2*pi; end
            if theta(2)-theta(1)<-pi, theta(2)=theta(2)+2*pi; end
            centerNear=centerAngle+2*pi*round((mean(theta)-centerAngle)/(2*pi));
            low=max(min(theta),centerNear-halfAngle);
            high=min(max(theta),centerNear+halfAngle);
            if high<=low, continue; end
            crossings=segmentCrossings(low-p.leakage.segmentStartAngleRad, ...
                high-p.leakage.segmentStartAngleRad,sectorBounds) ...
                +p.leakage.segmentStartAngleRad;
            cuts=[low,crossings,high];
            for piece=1:numel(cuts)-1
                fractions=sort(([cuts(piece),cuts(piece+1)]-theta(1))/diff(theta));
                midpoint=mean(fractions);
                len=norm(x(nodes(2),:)-x(nodes(1),:))*diff(fractions);
                sector=segmentIndex(mean(cuts(piece:piece+1))- ...
                    p.leakage.segmentStartAngleRad,sectorBounds);
                weight=len*[1-midpoint,midpoint];
                B.jointNodeWeight(jnt,local)=B.jointNodeWeight(jnt,local)+weight;
                B.jointSectorNodeWeight{sector}(jnt,local)= ...
                    B.jointSectorNodeWeight{sector}(jnt,local)+weight;
            end
        end
    end
    B.jointSupportLength=full(sum(B.jointNodeWeight,2));
    assert(all(B.jointSupportLength>0),'A joint band does not overlap the lining interface.');
else
    B.jointSupportLength=zeros(0,1);
end
[B.jointID,B.jointNodeIndex,jointWeight]=find(B.jointNodeWeight);
B.jointID=B.jointID(:); B.jointNodeIndex=B.jointNodeIndex(:);
jointWeight=jointWeight(:);
B.jointNormalizedWeight=jointWeight./B.jointSupportLength(B.jointID);
B.bodySectorNodeWeight=B.sectorNodeWeight;
B.jointSectorID=cell(p.leakage.segmentCount,1);
B.jointSectorNodeIndex=cell(p.leakage.segmentCount,1);
B.jointSectorNormalizedWeight=cell(p.leakage.segmentCount,1);
for sector=1:p.leakage.segmentCount
    sectorJointWeight=B.jointSectorNodeWeight{sector};
    B.bodySectorNodeWeight(:,sector)=B.bodySectorNodeWeight(:,sector)- ...
        full(sum(sectorJointWeight,1))';
    [jointID,nodeID,weight]=find(sectorJointWeight);
    jointID=jointID(:); nodeID=nodeID(:); weight=weight(:);
    B.jointSectorID{sector}=jointID;
    B.jointSectorNodeIndex{sector}=nodeID;
    B.jointSectorNormalizedWeight{sector}=weight./B.jointSupportLength(jointID);
end
assert(min(full(B.bodySectorNodeWeight(:)))>=-1e-8, ...
    'A joint band exceeds its physical segment support.');
B.bodySectorNodeWeight=max(B.bodySectorNodeWeight,0);
% A gasket band replaces the concrete body path over its physical width.
% Its integrated conductance is the total conductance of that band.
B.bodyW=B.drainW-full(sum(B.jointNodeWeight,1))';
assert(all(B.bodyW>=-1e-8),'Joint support exceeds the lining interface.');
B.bodyW=max(B.bodyW,0);
assert(max(abs(full(sum(B.bodySectorNodeWeight,2))-B.bodyW))<1e-8, ...
    'Body segments do not partition the remaining lining interface.');
B.far=[]; if p.hydroFar, B.far=G.map(unique([B.left;B.right;B.bottom])); end
B.topLocal=G.map(B.top);
B.hfar=x(G.ids(B.far),2)-G.zWT; % negative above WT, gravity balanced
B.Tdir=[];
if isinf(p.topHTC), B.Tdir=[B.Tdir;B.top]; end
if isinf(p.tunnelHTC), B.Tdir=[B.Tdir;find(B.tunnelW>0)]; end
B.Tdir=unique(B.Tdir);
if p.leakage.jointFreeze.enabled
    B.jointThermalPaths=build_joint_thermal_paths(S, ...
        p.leakage.jointAnglesRad,p.leakage.jointFreeze.quadraturePoints);
end
if isempty(p.leakage.jointAnglesRad)
    B.jointQuadrature=struct('phi',sparse(0,G.ns),'nj',0, ...
        'nsector',p.leakage.segmentCount,'jointID',zeros(0,1), ...
        'sectorID',zeros(0,1),'weight',zeros(0,1), ...
        'supportLength',zeros(0,1),'soilNodes',G.ids);
else
    B.jointQuadrature=build_joint_drain_quadrature(S,G.ids,p.leakage,B.segmentBoundaryAnglesRad);
end
B.soilDepth=G.xy(G.ids,2)-G.z0;
B.soilSelect=G.soilSelect;
end

function crossings=segmentCrossings(low,high,bounds)
% Local edge angles lie within about one circumference of the principal
% interval; copies on both sides handle the 0/2pi wrap without bin bias.
allBounds=reshape(bounds(:)+2*pi*(-1:1),[],1);
crossings=sort(allBounds(allBounds>low+1e-12 & allBounds<high-1e-12));
end

function sector=segmentIndex(angle,bounds)
sector=find(bounds<=mod(angle,2*pi),1,'last');
if isempty(sector), sector=numel(bounds); end
end

function w=lineWeights(x,n,dim)
[~,ix]=sort(x(n,dim)); n=n(ix); w=zeros(size(x,1),1);
if numel(n)>1, w=edgeWeights(x,[n(1:end-1),n(2:end)]); end
end

function w=edgeWeights(x,ed)
w=zeros(size(x,1),1);
if isempty(ed), return; end
L=sqrt(sum((x(ed(:,2),:)-x(ed(:,1),:)).^2,2));
w=accumarray(ed(:),repmat(L/2,2,1),[size(x,1) 1]);
end

function [X,t,info,G]=initialFields(G,B,p)
t=p.initial.startDay*p.day;
T=p.initial.surfaceT+p.initial.gradient*(G.xy(:,2)-G.z0);
h=G.xy(G.ids,2)-G.z0-waterTableDepthAt(t,p);
info=struct('mode',p.initial.mode, ...
    'note','Initial thermal and hydraulic transients are retained, not erased.', ...
    'priorRejectedSteps',0);
if strcmp(p.initial.mode,'checkpoint')
    C=load(p.initial.checkpointFile,'checkpoint'); c=C.checkpoint;
    assert(strcmp(c.version,'THM_v16_20260930_analytic_p1'), ...
        'Incompatible checkpoint version.');
    assert(isequal(c.pAll,G.xy) && isequal(c.soilNodes,G.ids) && isequal(c.tSoil,G.tri) ...
        && isequal(c.qsets,G.qsets) && isequal(c.ctri,G.ctri) && isequal(c.interface,G.interface),'Checkpoint mesh mismatch.');
    sameSignature=isequal(c.signature,signature(p));
    reintegrated=~sameSignature&&p.initial.allowNumericalReintegration ...
        &&thm_numerical_reintegration_compatible(c.signature,signature(p));
    assert(sameSignature||reintegrated,'THM:ConstitutiveMismatch', ...
        'Checkpoint constitutive parameters/closure differ. Start a new history for changed physics.');
    info.numericalReintegration=struct('used',reintegrated, ...
        'sourceSignature',c.signature,'targetSignature',signature(p), ...
        'scope','Same full state and installed reference; inventories re-integrated for a new conditional run.');
    if G.jointCycleSettings.enabled
        assert(isfield(c,'jointCycleMemory')&&isfield(c,'jointCycleAppliedSeverity'), ...
            'JointCycle:Checkpoint','Checkpoint lacks applied geometry and complete local cycle memory.');
        assert(c.jointCycleMemory.lastDay==c.time/p.day, ...
            'JointCycle:Checkpoint','Cycle memory is not from the checkpoint time.');
        G.jointCycleMemory=c.jointCycleMemory;
        validate_joint_cycle_memory(G.jointCycleMemory,8,G.jointThermalPaths.nq,G.jointCycleSettings);
        G=applyJointCycleGeometry(G,c.jointCycleAppliedSeverity);
    end
    T=c.X(1:G.N); h=c.X(G.N+(1:G.ns)); t=c.time;
    if isfield(c,'rejectedSteps')
        assert(isscalar(c.rejectedSteps) && c.rejectedSteps>=0 && ...
            c.rejectedSteps==fix(c.rejectedSteps),'Invalid checkpoint rejection count.');
        % A new run inherits the physical state, not the preceding run's
        % diagnostic counter. Only same-history recovery keeps that total.
        if p.initial.resumeHistory, info.priorRejectedSteps=c.rejectedSteps; end
    end
elseif strcmp(p.initial.mode,'fields')
    assert(numel(p.initial.T)==G.N && numel(p.initial.h)==G.ns,'Supply initial.T (all nodes) and initial.h (soil-node order).');
    T=p.initial.T(:); h=p.initial.h(:);
elseif ~isempty(p.initial.T) || ~isempty(p.initial.h)
    error('Use initial.mode=fields for supplied fields.');
end
assert(all(isfinite(T)) && all(T>-100 & T<100) && all(isfinite(h)));
X=[T;h];
if p.thm.enabled
    if strcmp(p.initial.mode,'checkpoint')
        assert(isfield(c,'mechanicalReference'),'THM checkpoint lacks installed reference.');
        G.referenceU=c.mechanicalReference.U;
        G.referenceT=c.mechanicalReference.T;
        X=c.X;
    else
        G.referenceT=T; G.referenceU=zeros(G.nm,1);
        U=thmInitialEquilibrium(T,h,G,p);
        G.referenceU=U; X=[X;U];
    end
    G.referenceStrainDiagnostics=soil_small_strain_diagnostics( ...
        G.strainOperator,G.referenceU(1:G.nuSoil));
    assert(G.referenceStrainDiagnostics.maxAbsPrincipalStrain<=p.thm.maximumPrincipalStrain ...
        &&G.referenceStrainDiagnostics.maxGradientNorm<=p.thm.maximumDisplacementGradient, ...
        'THM:StrainDomain','Installed soil equilibrium exceeds the small-strain domain.');
    info.mechanicalReference='Installed selfweight/hydrostatic/contact equilibrium; construction sequence not simulated.';
end
% No resetting measured fields to BC values here: changed forcing starts at first step.
info.initialFarBoundaryMismatch=0;
if ~isempty(B.far)
    boundaryHead=farHeadAt(t,G,B,p);
    info.initialFarBoundaryMismatch=max(abs(h(B.far)-boundaryHead));
end
end

function s=signature(p)
% Include every physical/boundary setting that affects the continuation.
s=thm_constitutive_signature(p);
end

function G=applyJointCycleGeometry(G,severity)
[G.jointLaw,info]=joint_cycle_channel_law(G.fixedJointLaw,severity,G.jointCycleSettings);
G.jointCycleAppliedSeverity=info.severity;
end

function c=makeCheckpoint(X,t,G,p,info,rejectedSteps)
c=struct('version','THM_v16_20260930_analytic_p1','X',X,'time',t,'pAll',G.xy,'soilNodes',G.ids, ...
    'tSoil',G.tri,'qsets',{G.qsets},'ctri',G.ctri,'interface',G.interface, ...
    'signature',signature(p),'parameters',p,'initialInfo',info, ...
    'rejectedSteps',rejectedSteps);
c.algorithmRevision='THM_v16_research_20261002_local_cycle_coarsening';
c.jointDamage=G.jointDamage;
if G.jointCycleSettings.enabled
    c.jointCycleMemory=G.jointCycleMemory;
    c.jointCycleAppliedSeverity=G.jointCycleAppliedSeverity;
end
if p.thm.enabled, c.mechanicalReference=struct('U',G.referenceU,'T',G.referenceT); end
end

function a=state(T,h,G,p,U,wantDerivatives)
% Diagnostic nodal phases; residual inventories are integrated inside SCVs.
if nargin<6, wantDerivatives=false; end
a=frozen_soil_constitutive(T(G.ids),h,p,wantDerivatives); a.T=T;
a.jointLaw=G.jointLaw;
Q=G.scv;
a.quadrature=frozen_soil_constitutive(Q.P*T(G.ids),Q.P*h,p,wantDerivatives);
if p.thm.enabled
    if nargin<5, U=G.referenceU; end
    a=thmState(a,T,h,U,G,p,wantDerivatives);
end
if G.soilConductivityField.enabled
    a=scale_soil_conductivity_state(a,G.soilConductivityField.nodalMultiplier);
    a.quadrature=scale_soil_conductivity_state( ...
        a.quadrature,G.soilConductivityField.quadratureMultiplier);
end
a.E=p.rhoC*p.cc*(G.Mconc*T);
a.E(G.ids)=a.E(G.ids)+Q.deposit*a.quadrature.energy;
a.backgroundEnergy=a.E;
a.W=Q.deposit*a.quadrature.mass;
a.elementK=Q.elementAverage*a.quadrature.K;
a.elementk=Q.elementAverage*a.quadrature.k;
if isfield(G,'jointThermalPaths') && p.leakage.jointStorage.enabled
    paths=G.jointThermalPaths; P=paths.interpolation;
    profile=reshape(P*T,paths.nq,paths.nj)';
    [channelE,channelDerivative,a.jointInventory]=joint_channel_inventory(profile, ...
        G.jointLaw.jointFreeze,G.jointLaw.jointStorage,G.jointLaw.jointConductance, ...
        p.liningThickness,p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
    a.E=a.E+P'*reshape(channelE',[],1);
    a.inventoryProfileDerivative=reshape(channelDerivative',[],1);
    a.W=a.W+G.channelWaterNodeMass;
end
end

function s=retention(h,p)
m=1-1/p.vgN; sr=p.thetaR/p.porosity;
s=sr+(1-sr)*(1+(p.vgAlpha*max(-h,0)).^p.vgN).^(-m);
end

function [R,D,a]=residual(X,old,t,dt,G,B,p,scaled,a0,a,wantDiagnostics)
T=X(1:G.N); h=X(G.N+(1:G.ns));
if nargin<10 || isempty(a), a=state(T,h,G,p,X(G.uIndex)); end
if nargin<11, wantDiagnostics=nargout>1; end
D=[];
a.waterTableDepth=waterTableDepthAt(t,p);
if nargin<9, a0=state(old(1:G.N),old(G.N+(1:G.ns)),G,p,old(G.uIndex)); end
ke=a.elementk; kt=sparse(G.ri,G.ci,G.kvals.*repmat(ke,9,1),G.N,G.N)+G.C;
trans=G.edgeElementPattern*a.elementK;
H=h-G.xy(G.ids,2); i=G.edges(:,1); j=G.edges(:,2); F=trans.*(H(i)-H(j));
qLeak=zeros(G.ns,1); qTop=qLeak;
channelHeatOn=p.advection && p.leakage.jointFreeze.enabled && ...
    p.leakage.jointFreeze.heatExchangeEnabled && p.tunnelDrain;
if wantDiagnostics
    leak=struct('body',0,'bySegment',zeros(p.leakage.segmentCount,1), ...
        'byJoint',zeros(numel(p.leakage.jointAnglesRad),1));
end
if p.tunnelDrain
    if wantDiagnostics || channelHeatOn, [qLeak,leak]=liningLeak(h,a,B,p);
    else, qLeak=liningLeak(h,a,B,p); end
elseif wantDiagnostics
    % Potential leakage uses this independently solved no-drain TH state.
    [~,leak]=liningLeak(h,a,B,p);
    leak.body=0; leak.bySegment(:)=0; leak.byJoint(:)=0;
end
if p.topSeepage
    qTop=p.topLeakConductance*B.topW(G.ids).*max(h,0);
end
rm=(a.W-a0.W)/dt+G.inc*F+qLeak+qTop;
qFar=zeros(G.ns,1); qFar(B.far)=-rm(B.far); % exact prescribed-head reaction, outward positive
re=(a.E-a0.E)/dt+kt*T;
if p.advection
    tw=T(G.ids);
    Fsub=soil_cvfe_flux(h,a.elementK,G.cvfe,'element');
    re=re+soil_cvfe_heat(T,Fsub,G.cvfe,p.rho,p.L+p.enthalpyReference_Jkg,p.cw);
    tfar=p.groundwaterInflowT+p.groundwaterInflowGradient*(G.xy(G.ids,2)-G.z0);
    Tadv=tw; Tadv(qFar<0)=tfar(qFar<0);
    re(G.ids)=re(G.ids)+p.rho*(p.L+p.enthalpyReference_Jkg+p.cw*tw).*(qLeak+qTop) ...
        +p.rho*(p.L+p.enthalpyReference_Jkg+p.cw*Tadv).*qFar;
end
channelExportCorrection=0;
if channelHeatOn
    [channelHeat,~,~,channelExportCorrection]=joint_channel_heat(T,G.soilSelect, ...
        B.jointThermalPaths,leak.nodeFlow,[],[],[],p.rho*p.cw);
    re=re+channelHeat;
end
[Ttop,Ttun]=thermalValues(t,p);
qHeat=-p.geothermalFlux*B.bottomW;
if isfinite(p.topHTC), qHeat=qHeat+p.topHTC*B.topW.*(T-Ttop); end
if isfinite(p.tunnelHTC), qHeat=qHeat+p.tunnelHTC*B.tunnelW.*(T-Ttun); end
re=re+qHeat;
qTdir=zeros(G.N,1); qTdir(B.Tdir)=-re(B.Tdir);
if wantDiagnostics
    D.water=sum(a.W); D.energy=sum(a.E);
    D.Qout=sum(qFar+qLeak+qTop); D.leak=sum(qLeak);
    D.leakBody=leak.body; D.leakBySegment=leak.bySegment;
    D.leakByJoint=leak.byJoint;
    D.jointFreezeRecord=[ones(1,numel(leak.byJoint)), ...
        nan(1,2*numel(leak.byJoint)),zeros(1,numel(leak.byJoint))];
    if isfield(leak,'freeze')
        f=leak.freeze;
        D.jointFreezeRecord=[f.factor(:)',f.minTemperature_C(:)', ...
            f.maxTemperature_C(:)',f.effectiveFrozenLengthFraction(:)'];
    end
    D.Eout=sum(qHeat+qTdir);
    if p.advection
        D.Eout=D.Eout+sum(p.rho*(p.L+p.enthalpyReference_Jkg+p.cw*tw).*(qLeak+qTop)+p.rho*(p.L+p.enthalpyReference_Jkg+p.cw*Tadv).*qFar);
    end
    D.Eout=D.Eout+channelExportCorrection;
    D.waterError=sum(a.W-a0.W)+dt*D.Qout;
    D.energyError=sum(a.E-a0.E)+dt*D.Eout;
    D.localWaterDefect=dt*(rm+qFar);
    D.localEnergyDefect=dt*(re+qTdir);
    D.minTrans=min(trans); D.waterReaction=qFar; D.F=F;
    D.jointHydraulicRecord=[leak.drivingHead(:)',leak.referenceHead(:)', ...
        leak.unrestrictedRate(:)',a.waterTableDepth,min(trans),sum(trans<0)];
end
R=dt*[re;rm];
if p.thm.enabled
    rmech=thmMechanicalResidual(T,h,X(G.uIndex),G,p,a);
    R=[R;rmech*p.thm.forceResidualScale];
    if wantDiagnostics
        D.maxVolumetricStrain=max(abs(a.elementVolumetricStrain));
        D.maxJointGap_m=max(a.contact.jointGap_m(:));
        D.maxJointPressure_Pa=max(a.contact.jointContactPressureGauss_Pa);
        D.mechanicalResidual_Nm=max(abs(rmech(G.mfree)));
        D.maxOmittedSlitVolumeFraction=a.maxOmittedSlitVolumeFraction;
        s=a.strainDiagnostics; sa=a.absoluteStrainDiagnostics;
        D.strainRecord=[t/p.day,s.maxAbsVolumetricStrain,s.maxAbsPrincipalStrain, ...
            s.maximumTensorShearStrain,s.maxGradientNorm,sa.maxAbsPrincipalStrain,sa.maxGradientNorm];
        c=a.contact; d=a.openingHydraulics;
        D.jointMechanicalRecord=[t/p.day,min(c.jointGap_m,[],1),mean(c.jointGap_m,1), ...
            max(c.jointGap_m,[],1),max(abs(c.jointSlip_m),[],1), ...
            max(c.jointContactPressure_Pa,[],1),c.boltForce_N_per_m', ...
            d.baseConductance',d.openConductance',a.effectiveJointConductance'];
    end
end
if scaled, R=R./G.rs; end
if p.thm.enabled
    R(G.N+G.ns+G.mfixed)=X(G.N+G.ns+G.mfixed)/p.thm.displacementScale_m;
end
% Boundary rows always dimensionless, same in numerical Jacobian and line search.
R(G.N+B.far)=(h(B.far)-farHeadAt(t,G,B,p))/10;
if isinf(p.topHTC), R(B.top)=T(B.top)-Ttop; end
if isinf(p.tunnelHTC), ids=find(B.tunnelW>0); R(ids)=T(ids)-Ttun; end
end

function [qLeak,ledger,QT,DH]=liningLeak(h,a,B,p)
% Consistent finite-width joint Robin quadrature; body path remains optional.
% The driving potential crosses the real lining thickness in gravity.
drivingHead=h-p.leakage.interiorHead+B.liningZDrop;
bodyFlux=p.liningK*B.bodyRadialFactor.*max(drivingHead,0);
body=B.bodyW.*bodyFlux;
if isfield(a,'waterTableDepth'), depth=a.waterTableDepth; else, depth=p.waterTableDepth; end
drainParameters=p;
drainParameters.leakage=a.jointLaw;
if p.thm.enabled
    drainParameters.leakage.jointConductance=a.effectiveJointConductance;
    drainParameters.leakage.jointFreeze.enabled=false;
end
[qJoint,ledger,Djoint,QT]=evaluate_joint_drain(h,a.T,B,drainParameters,depth);
if p.thm.enabled
    paths=B.jointThermalPaths;
    profile=reshape(paths.interpolation*a.T,paths.nq,paths.nj)';
    [~,~,ledger.freeze]=joint_freeze_conductance(profile,a.jointLaw.jointFreeze);
    if nargout>2
        unitParameters=drainParameters;
        unitParameters.leakage.jointConductance=ones(paths.nj,1);
        [~,unitLedger]=evaluate_joint_drain(h,a.T,B,unitParameters,depth);
        ledger.unitNodeFlow=unitLedger.nodeFlow;
        ledger.baseFlow=unitLedger.nodeFlow;
        ledger.dFactorT=a.jointConductanceT;
        QT=unitLedger.nodeFlow'*a.jointConductanceT;
    end
end
qLeak=body+qJoint;
ledger.body=sum(body);
ledger.bySegment=ledger.bySegment+B.bodySectorNodeWeight'*bodyFlux;
DH=Djoint+spdiags(p.liningK*B.bodyW.*B.bodyRadialFactor.* ...
    (drivingHead>0),0,numel(h),numel(h));
end

function [Ttop,Ttun]=thermalValues(t,p)
c=cos(2*pi*t/(p.yearDays*p.day));
Ttop=(p.topMax+p.topMin)/2-(p.topMax-p.topMin)/2*c;
Ttun=(p.tunnelMax+p.tunnelMin)/2-(p.tunnelMax-p.tunnelMin)/2*c;
phase=mod(t/p.day,p.yearDays);
if ~isempty(p.forcing.topAirC)
    Ttop=interp1(p.forcing.day,p.forcing.topAirC,phase,'linear');
end
if ~isempty(p.forcing.tunnelAirC)
    Ttun=interp1(p.forcing.day,p.forcing.tunnelAirC,phase,'linear');
end
end

function depth=waterTableDepthAt(t,p)
depth=p.waterTableDepth;
if ~isempty(p.forcing.waterTableDepth_m)
    phase=mod(t/p.day,p.yearDays);
    depth=interp1(p.forcing.day,p.forcing.waterTableDepth_m,phase,'linear');
end
end

function head=farHeadAt(t,G,B,p)
head=B.hfar;
if ~isempty(p.forcing.waterTableDepth_m)
    head=G.xy(G.ids(B.far),2)-G.z0-waterTableDepthAt(t,p);
end
end

function [cap,boundary]=calendarStepLimit(t,p)
% Numerical resolution only: resolve both sides of every scheduled change.
% The same twelve equal 365/12-day bins are used by leakage integration.
cap=p.dtMax; boundary=inf;
if isempty(p.solver.dtMaxByMonth_days), return; end
absoluteDay=t/p.day;
year=floor(absoluteDay/p.yearDays+1e-12);
edges=(0:12)*(p.yearDays/12);
phase=absoluteDay-year*p.yearDays;
month=find(phase<edges(2:end)-1e-9,1);
if isempty(month), month=1; year=year+1; end
cap=p.solver.dtMaxByMonth_days(month)*p.day;
boundary=(year*p.yearDays+edges(month+1))*p.day;
assert(boundary>t,'THM:CalendarStepCap','The next calendar boundary must advance time.');
end

function [X,info,cacheOut,aCurrent,a0]=step(old,t,dt,G,B,p,seed,cacheIn,a0)
% Adaptive modified Newton (chord Newton). A sparse analytic tangent is
% reused while it remains predictive.  Poor contraction forces a rebuild;
% a failed line search with a reused tangent gets one immediate fresh-tangent
% retry.  Therefore reuse changes cost, not the converged discrete equations.
X=old; J=[]; F=[]; Jage=inf; cacheOut=[]; R=[]; predictorEvaluations=0;
aCurrent=[];
if nargin>=8 && ~isempty(cacheIn) && p.solver.reuseAcrossSteps && ...
        abs(dt-cacheIn.dt)<1e-9*max(dt,cacheIn.dt)
    J=cacheIn.J; F=cacheIn.F; Jage=0;
end
if nargin<9 || isempty(a0)
    a0=state(old(1:G.N),old(G.N+(1:G.ns)),G,p,old(G.uIndex));
end
if p.solver.reuseConstitutiveState
    aCurrent=a0;
    % Old inventory must survive a coarsening event, but the CURRENT
    % residual and tangent must be re-evaluated with the new channel law.
    if ~isequal(a0.jointLaw,G.jointLaw),aCurrent=[];end
end
if nargin>=7 && ~isempty(seed) && thmAdmissible(seed,G,p)
    Rold=residual(old,old,t+dt,dt,G,B,p,true,a0,aCurrent,false);
    [Rseed,~,seedState]=residual(seed,old,t+dt,dt,G,B,p,true,a0,[],false);
    predictorEvaluations=2;
    if all(isfinite(Rseed)) && norm(Rseed,2)<norm(Rold,2)
        X=seed;
        aCurrent=seedState;
        if p.solver.reuseAcceptedResidual, R=Rseed; end
    elseif p.solver.reuseAcceptedResidual
        R=Rold;
    end
    if ~p.solver.reuseConstitutiveState, aCurrent=[]; end
end
info=struct('ok',false,'iter',0,'norm',inf,'reason','iteration limit', ...
    'jacobianBuilds',0,'linearSolves',0,'residualEvaluations',predictorEvaluations);
for k=1:p.solver.maxIter
    % A successful line search already evaluated this exact X/old/time/dt.
    % Retain its residual; recomputing it changes neither the tangent nor
    % acceptance and duplicates the expensive SCV constitutive evaluation.
    if isempty(R)
        [R,~,aCurrent]=residual(X,old,t+dt,dt,G,B,p,true,a0,aCurrent,false);
        info.residualEvaluations=info.residualEvaluations+1;
    end
    nr=norm(R,inf);
    info.iter=k; info.norm=nr;
    if ~isfinite(nr), info.reason='nonfinite residual'; return; end
    if nr<p.solver.tol
        info.ok=true; info.reason='converged';
        if ~isempty(F) && Jage<p.solver.jacobianReuse
            cacheOut=struct('J',J,'F',F,'dt',dt);
        end
        return
    end
    accepted=false;
    if isempty(J) || Jage>=p.solver.jacobianReuse
        [J,nfd]=jacobian(X,old,t+dt,dt,G,B,p,R,a0,aCurrent);
        F=factorize(J);
        info.jacobianBuilds=info.jacobianBuilds+1;
        info.residualEvaluations=info.residualEvaluations+nfd; Jage=0;
    end
    merit0=norm(R,2);
    for tangentAttempt=1:2
        dy=-(F\R); info.linearSolves=info.linearSolves+1; dx=G.xs.*dy;
        if any(~isfinite(dx)), info.reason='singular tangent'; return; end
        lam=min([1,p.solver.stepT/max(max(abs(dx(1:G.N))),eps), ...
            p.solver.stepH/max(max(abs(dx(G.N+(1:G.ns)))),eps)]);
        if p.thm.enabled, lam=min(lam,p.solver.stepU/max(max(abs(dx(G.uIndex))),eps)); end
        for ls=1:p.solver.lineSearch
            candidate=X+lam*dx;
            if all(candidate(1:G.N)>-100 & candidate(1:G.N)<100) && thmAdmissible(candidate,G,p)
                [rr,~,candidateState]=residual(candidate,old,t+dt,dt,G,B,p,true,a0,[],false);
                info.residualEvaluations=info.residualEvaluations+1;
                if norm(rr,2)<(1-1e-4*lam)*merit0
                    accepted=true; break
                end
            end
            lam=lam/2;
        end
        if accepted, break; end
        if Jage>0 && tangentAttempt==1
            % A stale chord must never be allowed to trigger a time-step cut.
            [J,nfd]=jacobian(X,old,t+dt,dt,G,B,p,R,a0,aCurrent);
            F=factorize(J);
            info.jacobianBuilds=info.jacobianBuilds+1;
            info.residualEvaluations=info.residualEvaluations+nfd; Jage=0;
        else
            break
        end
    end
    if ~accepted, info.reason='line search failed'; info.trialX=X; info.trialResidual=R; return; end
    contraction=norm(rr,inf)/max(nr,realmin);
    X=candidate; Jage=Jage+1;
    if p.solver.reuseConstitutiveState, aCurrent=candidateState; else, aCurrent=[]; end
    if p.solver.reuseAcceptedResidual, R=rr; else, R=[]; end
    if contraction>p.solver.jacobianRefreshRatio
        J=[]; F=[]; Jage=inf;
    end
end
info.trialX=X; info.trialResidual=R;
end

function F=factorize(J)
% Reuse the numeric LU factorization together with a reused chord tangent.
% The matrix fallback keeps compatibility with MATLAB releases lacking
% decomposition (backslash will then factorize for each solve).
if exist('decomposition','builtin') || exist('decomposition','file')
    F=decomposition(J,'lu');
else
    F=J;
end
end

function [J,nEval]=jacobian(X,old,t,dt,G,B,p,R,a0,aCurrent)
if nargin<9, a0=state(old(1:G.N),old(G.N+(1:G.ns)),G,p,old(G.uIndex)); end
if nargin<10,aCurrent=[];end
if ~p.thm.enabled || strcmp(p.solver.jacobian,'analytic')
    [J,nEval]=thJacobian(X,old,t,dt,G,B,p,R,a0,aCurrent); return
end
if nargin<9, a0=state(old(1:G.N),old(G.N+(1:G.ns)),G,p,old(G.uIndex)); end
v=zeros(numel(G.jrow),1); nEval=0;
for c=1:numel(G.groups)
    cols=G.groups{c}; [xp,delta]=thmPerturbation(X,cols,G,p);
    dr=residual(xp,old,t,dt,G,B,p,true,a0)-R; nEval=nEval+1;
    nz=G.jnz{c}; rows=G.jrow(nz); columns=G.jcol(nz);
    v(nz)=dr(rows)./delta(G.jloc{c}).*G.xs(columns);
end
J=sparse(G.jrow,G.jcol,v,numel(X),numel(X));
end

function [J,nEval]=thJacobian(X,old,t,dt,G,B,p,~,a0,aCurrent)
% Sparse tangent from analytic phase derivatives at positive SCV quadrature
% points, consistent concrete storage and exact physical subface operators.
if nargin<9, a0=state(old(1:G.N),old(G.N+(1:G.ns)),G,p,old(G.uIndex)); end
N=G.N; ns=G.ns; ne=size(G.tri,1); nf=size(G.edges,1);
T=X(1:N); h=X(N+(1:ns)); tw=T(G.ids);
if nargin<10 || isempty(aCurrent)
    a=state(T,h,G,p,X(G.uIndex),true);
else
    % The current residual state is exactly this X, mesh, material and
    % installed reference. Only SCV constitutive partials are missing.
    a=aCurrent;Tq=G.scv.P*T(G.ids);hq=G.scv.P*h;
    q=frozen_soil_constitutive(Tq,hq,p,true);
    if p.thm.enabled
        epsQ=G.qElementSelect*a.elementVolumetricStrain;
        q=thm_bulk_constitutive(Tq,hq,epsQ,p,q,true);
    end
    if G.soilConductivityField.enabled
        q=scale_soil_conductivity_state(q,G.soilConductivityField.quadratureMultiplier);
    end
    a.quadrature=q;
end
a.waterTableDepth=waterTableDepthAt(t,p);
Q=G.scv; aq=a.quadrature; nqp=size(Q.P,1);
% Exact sparse assembly of derivatives evaluated at the SAME SCV points.
if strcmp(p.solver.scvTangentAssembly,'element')
    lift=@(v) soil_scv_tangent_lift(v,Q.tangentAssembly);
    average=@(v) soil_scv_tangent_average(v,Q.tangentAssembly);
else
    lift=@(v) Q.deposit*spdiags(v,0,nqp,nqp)*Q.P;
    average=@(v) Q.elementAverage*spdiags(v,0,nqp,nqp)*Q.P;
end
dWT=lift(aq.dMassT); dWH=lift(aq.dMassH);
dETsoil=lift(aq.dEnergyT); dEHsoil=lift(aq.dEnergyH);
dkT=average(aq.dkT); dkH=average(aq.dkH);
dKT=average(aq.dKT); dKH=average(aq.dKH);
dWU=sparse(ns,G.nm); dEU=sparse(N,G.nm);
dkU=sparse(ne,G.nm); dKU=sparse(ne,G.nm);
if p.thm.enabled
    if strcmp(p.solver.scvTangentAssembly,'element')
        strainLift=@(v) soil_scv_tangent_strain(v,Q.tangentAssembly,G.msoil.V,'deposit');
        strainAverage=@(v) soil_scv_tangent_strain(v,Q.tangentAssembly,G.msoil.V,'average');
        dWU=[strainLift(aq.dMassEps),sparse(ns,G.nuRing)];
        dEU=G.soilSelect'*[strainLift(aq.dEnergyEps),sparse(ns,G.nuRing)];
        dkU=[strainAverage(aq.dkEps),sparse(ne,G.nuRing)];
        dKU=[strainAverage(aq.dKEps),sparse(ne,G.nuRing)];
    else
        strainMap=[G.poreDivergence,sparse(nqp,G.nuRing)];
        dWU=Q.deposit*spdiags(aq.dMassEps,0,nqp,nqp)*strainMap;
        dEU=G.soilSelect'*Q.deposit*spdiags(aq.dEnergyEps,0,nqp,nqp)*strainMap;
        dkU=Q.elementAverage*spdiags(aq.dkEps,0,nqp,nqp)*strainMap;
        dKU=Q.elementAverage*spdiags(aq.dKEps,0,nqp,nqp)*strainMap;
    end
end
ke=a.elementk;
kt=sparse(G.ri,G.ci,G.kvals.*repmat(ke,9,1),N,N)+G.C;
base=reshape(G.kvals,ne,9); localAction=zeros(ne,3);
for i=1:3
    for j=1:3
        localAction(:,i)=localAction(:,i)+base(:,(i-1)*3+j).*T(G.tri(:,j));
    end
end
action=sparse(G.tri(:),repmat((1:ne)',3,1),localAction(:),N,ne);
thermalT=action*dkT*G.soilSelect;
thermalH=action*dkH;
thermalU=action*dkU;
H=h-G.xy(G.ids,2); i=G.edges(:,1); j=G.edges(:,2);
headDiff=H(i)-H(j); trans=G.edgeElementPattern*a.elementK;
edgeProperty=spdiags(headDiff,0,nf,nf)*G.edgeElementPattern;
dFT=edgeProperty*dKT;
dFU=edgeProperty*dKU;
dFH=edgeProperty*dKH + ...
    sparse((1:nf)',i,trans,nf,ns)-sparse((1:nf)',j,trans,nf,ns);
F=trans.*headDiff;
% Linearized head-dependent lining and surface drainage.
qLeak=zeros(ns,1); QT=sparse(ns,N); Dq=sparse(ns,ns); QU=sparse(ns,G.nm);
if p.tunnelDrain
    [qLeak,leak,QT,Dq]=liningLeak(h,a,B,p);
    if p.thm.enabled
        QU=[sparse(ns,G.nuSoil),leak.unitNodeFlow'*a.jointConductanceU];
    end
end
qTop=zeros(ns,1);
if p.topSeepage
    qTop=p.topLeakConductance*B.topW(G.ids).*max(h,0);
    Dq=Dq+spdiags(p.topLeakConductance*B.topW(G.ids).*(h>0),0,ns,ns);
end
dRmTlocal=dWT+dt*(G.inc*dFT);
dRmT=dRmTlocal*G.soilSelect+dt*QT;
dRmH=dWH+dt*(G.inc*dFH+Dq);
dRmU=dWU+dt*(G.inc*dFU+QU);
dReT=p.rhoC*p.cc*G.Mconc+G.soilSelect'*dETsoil*G.soilSelect+dt*(kt+thermalT);
if isfield(a,'inventoryProfileDerivative')
    P=G.jointThermalPaths.interpolation;
    dReT=dReT+P'*spdiags(a.inventoryProfileDerivative,0,size(P,1),size(P,1))*P;
end
dReH=G.soilSelect'*dEHsoil+dt*thermalH;
dReU=dEU+dt*thermalU;
if p.advection
    [Fsub,DsubH,DsubK]=soil_cvfe_flux(h,a.elementK,G.cvfe,'element');
    [~,heatT,heatF]=soil_cvfe_heat(T,Fsub,G.cvfe,p.rho,p.L+p.enthalpyReference_Jkg,p.cw);
    subT=DsubK*dKT*G.soilSelect;
    subH=DsubH+DsubK*dKH;
    subU=DsubK*dKU;
    dReT=dReT+dt*(heatT+heatF*subT);
    dReH=dReH+dt*(heatF*subH);
    dReU=dReU+dt*(heatF*subU);
    rm=(a.W-a0.W)/dt+G.inc*F+qLeak+qTop;
    qFar=zeros(ns,1); qFar(B.far)=-rm(B.far);
    tfar=p.groundwaterInflowT+p.groundwaterInflowGradient*(G.xy(G.ids,2)-G.z0);
    Tadv=tw; Tadv(qFar<0)=tfar(qFar<0);
    wFar=p.rho*(p.L+p.enthalpyReference_Jkg+p.cw*Tadv);
    wDrain=p.rho*(p.L+p.enthalpyReference_Jkg+p.cw*tw);
    farMask=spdiags(double(ismember((1:ns)',B.far)),0,ns,ns);
    dBoundT=dt*spdiags(p.rho*p.cw*(qLeak+qTop+qFar.*(qFar>=0)),0,ns,ns)*G.soilSelect ...
        +dt*spdiags(wDrain,0,ns,ns)*QT ...
        -spdiags(wFar,0,ns,ns)*farMask*dRmT;
    dBoundH=dt*spdiags(wDrain,0,ns,ns)*Dq ...
        -spdiags(wFar,0,ns,ns)*farMask*dRmH;
    dBoundU=dt*spdiags(wDrain,0,ns,ns)*QU ...
        -spdiags(wFar,0,ns,ns)*farMask*dRmU;
    dReT=dReT+G.soilSelect'*dBoundT;
    dReH=dReH+G.soilSelect'*dBoundH;
    dReU=dReU+G.soilSelect'*dBoundU;
end
heatSlope=zeros(N,1);
if p.advection && p.tunnelDrain && p.leakage.jointFreeze.enabled && ...
        p.leakage.jointFreeze.heatExchangeEnabled
    [~,channelT,channelH]=joint_channel_heat(T,G.soilSelect,B.jointThermalPaths, ...
        leak.nodeFlow,leak.baseFlow,leak.headSlope,leak.dFactorT,p.rho*p.cw, ...
        leak.entryHeatHeadDerivative);
    dReT=dReT+dt*channelT; dReH=dReH+dt*channelH;
    if p.thm.enabled
        paths=B.jointThermalPaths; nj=paths.nj; nq=paths.nq; np=nj*nq;
        P=paths.interpolation;
        first=sparse((1:nq:np)',(1:nj)',1,np,nj);
        repeat=kron(speye(nj),ones(nq,1));
        difference=kron(speye(nj),-speye(nq)+sparse((2:nq)',(1:nq-1)',1,nq,nq));
        unit=leak.unitNodeFlow;
        actionU=first*spdiags(unit*tw,0,nj,nj) ...
            +spdiags(difference*(P*T),0,np,np)*repeat*spdiags(full(sum(unit,2)),0,nj,nj);
        channelU=-p.rho*p.cw*P'*actionU*a.jointConductanceU;
        dReU=dReU+dt*[sparse(N,G.nuSoil),channelU];
    end
end
if isfinite(p.topHTC), heatSlope=heatSlope+p.topHTC*B.topW; end
if isfinite(p.tunnelHTC), heatSlope=heatSlope+p.tunnelHTC*B.tunnelW; end
dReT=dReT+dt*spdiags(heatSlope,0,N,N);
J=[dReT,dReH;dRmT,dRmH];
if p.thm.enabled
    % The mechanical tangent uses the contact state already evaluated above.
    Jmm=[G.mKss,G.mKrs';G.mKrs,a.contact.ringTangent_N_per_m2];
    Tq=Q.P*tw; hq=Q.P*h; pw=p.rho*p.g*hq;
    pi=pw-p.rho*p.L/p.Tm*min(Tq,0);
    if strcmp(p.thm.phasePressure,'bishopIce')
        pressureT=aq.slT.*pw+aq.siT.*pi-aq.si*p.rho*p.L/p.Tm.*(Tq<0);
        pressureH=aq.slH.*pw+aq.siH.*pi+p.rho*p.g*(aq.sl+aq.si);
    else
        pressureT=zeros(nqp,1); pressureH=p.rho*p.g*ones(nqp,1);
    end
    soilT=-G.msoil.A_temperature*G.msoil.elementMeanNodal*G.soilSelect ...
        -G.msoil.A_pressure*average(pressureT)*G.soilSelect;
    soilH=-G.msoil.A_pressure*average(pressureH);
    JmT=[soilT;-G.ring.thermalLoadMap];
    JmH=[soilH;sparse(G.nuRing,ns)];
    J=[J,[dReU;dRmU];p.thm.forceResidualScale*[JmT,JmH,Jmm]];
end
nx=size(J,1);
J=spdiags(1./G.rs,0,nx,nx)*J*spdiags(G.xs,0,nx,nx);
if p.thm.enabled
    rows=N+ns+G.mfixed; J(rows,:)=0; J(rows,rows)=speye(numel(rows));
end
% Dirichlet rows overwrite the physical residual and must overwrite its
% tangent as well. The hydraulic scaled unknown is 10 h, yielding unity.
if ~isempty(B.far)
    rows=N+B.far;
    J(rows,:)=0; J(rows,rows)=speye(numel(rows));
end
if isinf(p.topHTC), J(B.top,:)=0; J(B.top,B.top)=speye(numel(B.top)); end
if isinf(p.tunnelHTC)
    rows=find(B.tunnelW>0);
    J(rows,:)=0; J(rows,rows)=speye(numel(rows));
end
nEval=0;
end

function q=darcy(X,G,p)
a=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex)); h=X(G.N+(1:G.ns));
gradH=zeros(numel(G.area),2);
for j=1:3, gradH=gradH+G.grad(:,:,j).*h(G.local(:,j)); end
gradH(:,2)=gradH(:,2)-1;
q.elemFlux=-a.elementK.*gradH;
q.elemCentroid=(G.xy(G.tri(:,1),:)+G.xy(G.tri(:,2),:)+G.xy(G.tri(:,3),:))/3;
q.head=h-G.xy(G.ids,2); [q.subfaceFlux_m2s,~,~,q.cvfe]=soil_cvfe_flux(h,a.elementK,G.cvfe,'element');
q.note='Element Darcy velocity and actual median-dual subface flux; phase-integrated conductivity.';
end

function m=mechanics(X,X0,G,~,p)
% Auxiliary elastic equilibrium about initial state. No feedback to TH.
% Correct 3D isotropic eigenstrain under plane strain: eigenstress=3*Kbulk*eps0.
n=G.N; K=sparse(2*n,2*n); f=zeros(2*n,1); f0=f;
a=state(X(1:n),X(n+1:end),G,p); a0=state(X0(1:n),X0(n+1:end),G,p);
D=elasticD(p.mechanics.Esoil,p.mechanics.nuSoil);
sigma0=p.mechanics.Esoil/(1-2*p.mechanics.nuSoil)*[1;1;0];
for e=1:size(G.tri,1)
    ids=G.tri(e,:); loc=G.local(e,:); B=bmatrix(squeeze(G.grad(e,:,:))'); dof=udofs(ids);
    K(dof,dof)=K(dof,dof)+B'*D*B*G.area(e);
    et=p.mechanics.alphaSoil*mean(X(ids)); et0=p.mechanics.alphaSoil*mean(X0(ids));
    % Optional empirical extra strain. Off by default; not an ice-lens model.
    et=et+p.mechanics.frostEfficiency*(1000/917-1)*mean(a.thetaI(loc))/3;
    et0=et0+p.mechanics.frostEfficiency*(1000/917-1)*mean(a0.thetaI(loc))/3;
    pp=p.mechanics.biot*mean(a.sl(loc).*X(n+loc))*p.rho*p.g;
    pp0=p.mechanics.biot*mean(a0.sl(loc).*X0(n+loc))*p.rho*p.g;
    f(dof)=f(dof)+B'*(sigma0*et+[1;1;0]*pp)*G.area(e);
    f0(dof)=f0(dof)+B'*(sigma0*et0+[1;1;0]*pp0)*G.area(e);
end
D=elasticD(p.mechanics.Econc,p.mechanics.nuConc);
sigma0=p.mechanics.Econc/(1-2*p.mechanics.nuConc)*p.mechanics.alphaConc*[1;1;0];
for set=1:numel(G.qsets)
    qs=G.qsets{set};
    for e=1:size(qs,1)
        ids=qs(e,:); dof=udofs(ids); x=G.xy(ids,:);
        for xi=[-1 1]/sqrt(3)
            for eta=[-1 1]/sqrt(3)
                [N,d]=qshape(xi,eta); J=d*x; dn=(J\d)'; B=bmatrix(dn); da=abs(det(J));
                K(dof,dof)=K(dof,dof)+B'*D*B*da;
                f(dof)=f(dof)+B'*sigma0*(N*X(ids))*da;
                f0(dof)=f0(dof)+B'*sigma0*(N*X0(ids))*da;
            end
        end
    end
end
if ~isempty(G.ctri)
    [ar,gr]=triGeometry(G.xy,G.ctri);
    for e=1:size(G.ctri,1)
        ids=G.ctri(e,:); dof=udofs(ids); B=bmatrix(squeeze(gr(e,:,:))');
        K(dof,dof)=K(dof,dof)+B'*D*B*ar(e);
        f(dof)=f(dof)+B'*sigma0*mean(X(ids))*ar(e);
        f0(dof)=f0(dof)+B'*sigma0*mean(X0(ids))*ar(e);
    end
end
for e=1:size(G.interface,1)
    ids=G.interface(e,:); x=G.xy(ids,:); tangent=(x(2,:)-x(1,:))'; L=norm(tangent); tangent=tangent/L;
    normal=[tangent(2);-tangent(1)]; Dk=p.mechanics.kn*(normal*normal')+p.mechanics.kt*(tangent*tangent');
    dof=udofs(ids); Me=zeros(8);
    for xi=[-1 1]/sqrt(3)
        N=.5*[1-xi,1+xi]; A=[-N(1)*eye(2),-N(2)*eye(2),N(1)*eye(2),N(2)*eye(2)];
        Me=Me+A'*Dk*A*L/2;
    end
    K(dof,dof)=K(dof,dof)+Me;
end
x=G.xy; tol=1e-8;
side=find(abs(x(:,1)-min(x(:,1)))<tol | abs(x(:,1)-max(x(:,1)))<tol);
bottom=find(abs(x(:,2)-max(x(:,2)))<tol);
fixed=unique([2*side-1;2*bottom]); free=setdiff((1:2*n)',fixed);
u0=zeros(2*n,1); u=u0; u0(free)=K(free,free)\f0(free); u(free)=K(free,free)\f(free);
m.U=u-u0; m.UinitialEquilibrium=u0; m.UabsoluteReference=u;
m.note='One-way illustrative elastic response, relative to initial equilibrium; no geostatic/excavation or ice-lens constitutive model.';
m.freeResidual=norm(K(free,:)*m.U-(f(free)-f0(free)),inf)/max(1,norm(f(free)-f0(free),inf));
end

function D=elasticD(E,nu)
D=E/((1+nu)*(1-2*nu))*[1-nu nu 0;nu 1-nu 0;0 0 (1-2*nu)/2];
end

function B=bmatrix(d)
B=zeros(3,2*size(d,1));
for i=1:size(d,1), B(:,2*i-1:2*i)=[d(i,1) 0;0 d(i,2);d(i,2) d(i,1)]; end
end

function dof=udofs(ids)
dof=reshape([2*ids(:)'-1;2*ids(:)'],1,[]);
end

function plotOutput(out,G)
vals={out.T,out.Psoil/1000,out.state.thetaW,out.state.thetaI};
names={'Temperature (C)','Liquid pore pressure (kPa)','Liquid water content','Ice content (equal-density TH)'};
figure('Name','Conservative TH v5','Color','w');
for k=1:4
    subplot(2,2,k); v=nan(G.N,1);
    if k==1, v=vals{k}; else, v(G.ids)=vals{k}; end
    patch('Faces',G.tri,'Vertices',G.xy,'FaceVertexCData',v,'FaceColor','interp','EdgeColor','none'); hold on;
    if k==1
        for s=1:numel(G.qsets), patch('Faces',G.qsets{s},'Vertices',G.xy,'FaceVertexCData',v,'FaceColor','interp','EdgeColor','none'); end
        if ~isempty(G.ctri), patch('Faces',G.ctri,'Vertices',G.xy,'FaceVertexCData',v,'FaceColor','interp','EdgeColor','none'); end
    end
    axis equal tight; set(gca,'YDir','reverse'); colorbar; xlabel('Y (m)'); ylabel('Z down (m)'); title(names{k});
end
figure('Name','TH conservation','Color','w');
subplot(2,1,1); plot(out.history(:,1),cumsum(out.history(:,8))); ylabel('Cumulative water defect (m2)'); grid on;
subplot(2,1,2); plot(out.history(:,1),cumsum(out.history(:,9))); ylabel('Cumulative energy defect (J/m)'); xlabel('Day'); grid on;
end

function report=selftest()
% Small independent-invariant tests; results returned and printed.
p=defaults(); p.thm.enabled=false; p.plot=false; p.saveOutput=false; p.verbose=false;
p.mechanics.enable=false; p.hydroFar=false; p.topSeepage=false; p.tunnelDrain=false;
p.topHTC=15; p.geothermalFlux=0; p.topMin=2; p.topMax=2;
S=struct('pAll',[0 0;1 0;0 1;1 1],'tSoil',[1 2 3;2 4 3]);
G=geometry(S,p); B=boundaries(G,S,p);
% 1 affine gradient and nodal volumes
assert(norm(G.vol-[1;2;2;1]/6,inf)<1e-14);
h=G.xy(:,1)+2*G.xy(:,2); qg=zeros(2,2);
for j=1:3, qg=qg+G.grad(:,:,j).*h(G.local(:,j)); end
assert(norm(qg-repmat([1 2],2,1),inf)<1e-13);
% 2 phase closure, monotone mass vs h, and energy vs T
for T=[-10 -1 -.01 0 5]
    hh=linspace(-100,10,4)'; a=state(T*ones(4,1),hh,G,p);
    assert(all(a.thetaW>=0 & a.thetaI>=0 & a.thetaA>=-1e-14));
    assert(max(abs(a.thetaW+a.thetaI+a.thetaA-p.porosity))<1e-13);
    assert(max(abs(a.sl./(1-a.si)-retention(hh,p)))<1e-12);
    ap=state(T*ones(4,1),hh+1e-4,G,p); assert(all(ap.W>=a.W));
end
% 3 gravity-balanced hydrostatic field, no flux, closed isothermal equilibrium
p.topMin=2; p.topMax=2; X=[2*ones(4,1);G.xy(:,2)-5];
[R,D]=residual(X,X,86400,86400,G,B,p,true);
assert(norm(R,inf)<1e-12); assert(max(abs(D.F))<1e-15);
% 4 colored Jacobian against uncompressed perturbations in mixed frozen state
Xj=[-.04;.15;-.02;.05;-3;-2;-2;-1]; old=Xj;
R=residual(Xj,old,100,100,G,B,p,true); J=jacobian(Xj,old,100,100,G,B,p,R);
Jd=zeros(numel(Xj),numel(Xj));
for j=1:numel(Xj)
    d=p.solver.fd*max(abs(Xj(j)),G.xs(j)); z=Xj; z(j)=z(j)+d;
    Jd(:,j)=(residual(z,old,100,100,G,B,p,true)-R)/d*G.xs(j);
end
jacerr=norm(full(J)-Jd,inf)/max(1,norm(Jd,inf)); assert(jacerr<1e-7);
% 5 cooling a sealed specimen: nonzero ice, conserved water and heat budget
p.topMin=-2; p.topMax=-2; X=[.2*ones(4,1);G.xy(:,2)-5]; astart=state(X(1:4),X(5:8),G,p); W0=sum(astart.W);
t=0; maxW=0; maxE=0; steps=80;
for k=1:steps
    [z,info]=step(X,t,3600,G,B,p); assert(info.ok,'Sealed freezing nonlinear solve failed.');
    [~,d]=residual(z,X,t+3600,3600,G,B,p,false);
    maxW=max(maxW,abs(d.waterError)); maxE=max(maxE,abs(d.energyError));
    X=z; t=t+3600;
end
a=state(X(1:4),X(5:8),G,p); assert(max(a.thetaI)>.001);
assert(abs(sum(a.W)-W0)<1e-7); assert(maxW<1e-7 && maxE<10);
% 6 open-system advection/reaction accounting and Dirichlet heat accounting
p.hydroFar=true; p.topHTC=inf; p.topMin=1; p.topMax=1; B=boundaries(G,S,p);
B.hfar=B.hfar+.5*G.xy(G.ids(B.far),1);
z=[ones(4,1);G.xy(:,2)-5+.5*G.xy(:,1)]; [znew,info]=step(z,0,10000,G,B,p); assert(info.ok);
[~,d]=residual(znew,z,10000,10000,G,B,p,false);
assert(sum(abs(d.waterReaction))>1e-12);
assert(abs(d.waterError)<1e-7 && abs(d.energyError)<10);
% 7 mechanical reference equilibrium and corrected plane-strain thermal patch
p.mechanics.enable=true; m=mechanics(znew,z,G,S,p); assert(m.freeResidual<1e-10);
nu=.35; DD=elasticD(80e6,nu); q=(1+nu)*DD*[1;1;0];
assert(norm(q-80e6/(1-2*nu)*[1;1;0],inf)<1e-6);
% 8 real-interface leakage partition and finite-width gasket conductance
root=fileparts(mfilename('fullpath'));
mesh=load(fullfile(root,'example_mesh_verified.mat'));
p=defaults(); p.thm.enabled=false; p.mechanics.enable=false;
GG=geometry(mesh,p); BB=boundaries(GG,mesh,p);
segmentLengths=full(sum(BB.sectorNodeWeight,1));
assert(max(abs(segmentLengths-mean(segmentLengths)))<1e-8);
hh=20*ones(GG.ns,1)-BB.liningZDrop;
aa=state(2*ones(GG.N,1),hh,GG,p);
[leakNode,ledger]=liningLeak(hh,aa,BB,p);
bodyExpected=sum(BB.drainW.*BB.bodyRadialFactor)*p.liningK*20;
assert(abs(sum(leakNode)-bodyExpected)<1e-12);
assert(abs(sum(ledger.bySegment)-bodyExpected)<1e-12);
assert(abs(mean(BB.bodyRadialFactor(BB.drainW>0))- ...
    1/(6*log(6/5.55)))<1e-10,'Annular radial conductance is incorrect.');
balancedHead=-BB.liningZDrop;
balancedLeak=liningLeak(balancedHead,aa,BB,p);
assert(sum(balancedLeak)<1e-14,'Equal total head must cause zero lining leakage.');
interfaceNodes=unique(mesh.soilTunnelCoupleTable(:,[2 3]));
[~,crown]=min(mesh.pAll(interfaceNodes,2));
[~,invert]=max(mesh.pAll(interfaceNodes,2));
crownLocal=GG.map(interfaceNodes(crown));
invertLocal=GG.map(interfaceNodes(invert));
assert(BB.liningZDrop(crownLocal)>0 && BB.liningZDrop(invertLocal)<0);
gravityOnly=liningLeak(zeros(GG.ns,1),aa,BB,p);
assert(gravityOnly(crownLocal)>0 && gravityOnly(invertLocal)==0);
p.liningK=0; p.leakage.jointAnglesRad=0;
p.leakage.jointConductance=2e-8;
p.leakage.jointEntryHead=1; p.leakage.interiorHead=.5;
hh=3.5*ones(GG.ns,1);
for width=[.02 .5 1.5]
    p.leakage.jointWidth_m=width;
    BB=boundaries(GG,mesh,p);
    [leakNode,ledger]=liningLeak(hh,aa,BB,p);
    assert(abs(sum(leakNode)-4e-8)<1e-14);
    assert(abs(ledger.byJoint-4e-8)<1e-14);
    assert(abs(sum(ledger.bySegment)-4e-8)<1e-14);
    assert(abs(sum(BB.bodyW)+sum(BB.jointSupportLength)-sum(BB.drainW))<1e-9);
    assert(ledger.bySegment(1)>0 && ledger.bySegment(end)>0, ...
        'A finite-width boundary joint must count toward both adjacent segments.');
    assert(abs(ledger.bySegment(1)+ledger.bySegment(end)-ledger.byJoint)<1e-14);
end
p.liningK=1e-10;
[mixedLeak,mixedLedger]=liningLeak(hh,aa,BB,p);
assert(abs(sum(mixedLeak)-mixedLedger.body-sum(mixedLedger.byJoint))<1e-12);
assert(abs(sum(mixedLedger.bySegment)-sum(mixedLeak))<1e-12, ...
    'Body and finite-width joint rates must partition all segment rates.');
p.liningK=0;
hh=ones(GG.ns,1); leakNode=liningLeak(hh,aa,BB,p);
assert(all(leakNode==0));
p.leakage.jointEntryHead=0; p.leakage.interiorHead=0;
p.leakage.jointWidth_m=.02; p.leakage.jointAnglesRad=3*pi/2;
BB=boundaries(GG,mesh,p);
[~,crownJoint]=liningLeak(zeros(GG.ns,1),aa,BB,p);
assert(crownJoint.byJoint>0,'Gravity must drive a crown joint.');
p.leakage.jointAnglesRad=pi/2;
BB=boundaries(GG,mesh,p);
[~,invertJoint]=liningLeak(zeros(GG.ns,1),aa,BB,p);
assert(invertJoint.byJoint==0,'Gravity must oppose an invert joint.');
% 9 timestep regression: finite-difference Jacobian cost cannot trap
% successful Newton steps at seconds-long increments.
p=defaults(); p.thm.enabled=false; p.solver.maxIter=60;
effortA=struct('iter',12,'residualEvaluations',50);
effortB=struct('iter',12,'residualEvaluations',1000);
dtA=acceptedTimeStep(.02*p.day,.02*p.day,effortA,p);
dtB=acceptedTimeStep(.02*p.day,.02*p.day,effortB,p);
assert(dtA==dtB && dtA>.02*p.day);
moderate=struct('iter',34,'residualEvaluations',1000);
recovered=.002*p.day;
for k=1:100
    recovered=acceptedTimeStep(recovered,recovered,moderate,p);
end
assert(recovered>.1*p.day,'Accepted moderate steps did not recover.');
hard=struct('iter',50,'residualEvaluations',50);
assert(acceptedTimeStep(.5*p.day,.5*p.day,hard,p)<.5*p.day);
assert(acceptedTimeStep(.5*p.day,.01*p.day,hard,p)>=.5*p.day, ...
    'An output-boundary remainder must not shrink the proposed step.');
% 10 real-mesh tangent: simultaneous freezing/thawing, wet/dry upper
% seepage, open/closed finite-width joints, and far-field reactions.
p=defaults(); p.thm.enabled=false; p.plot=false; p.saveOutput=false; p.verbose=false;
p.mechanics.enable=false;
p.leakage.jointAnglesRad=[0;pi/2];
p.leakage.jointWidth_m=[.25;.25];
p.leakage.jointConductance=[2e-8;1e-8];
p.leakage.jointEntryHead=[0;100];
mesh=load(fullfile(root,'example_mesh_repaired_interface.mat'));
GG=geometry(mesh,p); BB=boundaries(GG,mesh,p);
xy=GG.xy; z=xy(:,2)-GG.z0;
TT=-2+.11*z+.3*sin(xy(:,1)/7);
hh=z(GG.ids)-5+2*sin(xy(GG.ids,1)/6);
topLocal=GG.map(BB.top);
hh(topLocal)=.3+.8*sin(xy(BB.top,1)/3);
XX=[TT;hh]; XXold=XX;
XXold(GG.N+BB.far)=hh(BB.far)+10*sin(xy(GG.ids(BB.far),2)/3);
dt=3600; time=19*365*p.day+dt;
[RR,DD]=residual(XX,XXold,time,dt,GG,BB,p,true);
assert(any(TT(GG.ids)<0) && any(TT(GG.ids)>0));
assert(any(hh(topLocal)>0) && any(hh(topLocal)<0));
assert(any(DD.waterReaction(BB.far)>0) && any(DD.waterReaction(BB.far)<0));
assert(any(DD.leakByJoint>0) && any(DD.leakByJoint==0));
JJ=jacobian(XX,XXold,time,dt,GG,BB,p,RR);
v=zeros(numel(GG.jrow),1);
for c=1:numel(GG.groups)
    cols=GG.groups{c}; delta=p.solver.fd*max(abs(XX(cols)),GG.xs(cols));
    xp=XX; xp(cols)=xp(cols)+delta;
    dr=residual(xp,XXold,time,dt,GG,BB,p,true)-RR;
    nz=GG.jnz{c}; ii=GG.jrow(nz); jj=GG.jcol(nz);
    v(nz)=dr(ii)./delta(GG.jloc{c}).*GG.xs(jj);
end
JJfd=sparse(GG.jrow,GG.jcol,v,numel(XX),numel(XX));
realMeshJacerr=norm(JJ-JJfd,inf)/max(1,norm(JJfd,inf));
assert(realMeshJacerr<1e-5,'Real-mesh sparse tangent differs from coloured FD.');
% 11 Frozen-joint nonlocal tangent. Temperature DOFs inside the lining
% must influence remote soil mass, drainage enthalpy, and far reactions.
p.leakage.jointFreeze.enabled=true;
p.leakage.jointFreeze.flowWeightedRadius_m=1e-6;
BB=boundaries(GG,mesh,p); GG=jointCouplingPattern(GG,BB);
GG.jointThermalPaths=BB.jointThermalPaths;
Q=BB.jointQuadrature; nq=BB.jointThermalPaths.nq; nj=BB.jointThermalPaths.nj;
Jmap=sparse((1:numel(Q.weight))',Q.jointID,1,numel(Q.weight),nj);
[~,~,inventory]=joint_channel_inventory(zeros(nj,nq),p.leakage.jointFreeze, ...
    p.leakage.jointStorage,p.leakage.jointConductance,p.liningThickness, ...
    p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
GG.channelWaterNodeMass=Q.phi'*spdiags(Q.weight,0,numel(Q.weight), ...
    numel(Q.weight))*Jmap*inventory.waterVolume_m2;
TT(GG.vconc>0)=-.05;
XX=[TT;hh]; XXold=XX;
XXold(GG.N+BB.far)=hh(BB.far)+10*sin(xy(GG.ids(BB.far),2)/3);
[RR,~]=residual(XX,XXold,time,dt,GG,BB,p,true);
JJ=jacobian(XX,XXold,time,dt,GG,BB,p,RR);
v=zeros(numel(GG.jrow),1);
for c=1:numel(GG.groups)
    cols=GG.groups{c}; delta=p.solver.fd*max(abs(XX(cols)),GG.xs(cols));
    xp=XX; xp(cols)=xp(cols)+delta;
    dr=residual(xp,XXold,time,dt,GG,BB,p,true)-RR;
    nz=GG.jnz{c}; ii=GG.jrow(nz); jj=GG.jcol(nz);
    v(nz)=dr(ii)./delta(GG.jloc{c}).*GG.xs(jj);
end
JJfd=sparse(GG.jrow,GG.jcol,v,numel(XX),numel(XX));
freezeJacerr=norm(JJ-JJfd,inf)/max(1,norm(JJfd,inf));
if freezeJacerr>=1e-5
    [erow,ecol,evalue]=find(JJ-JJfd); [~,worst]=max(abs(evalue));
    fprintf('Frozen tangent error %.6g; worst row %d col %d: analytic %.8g, FD %.8g.\n', ...
        freezeJacerr,erow(worst),ecol(worst),full(JJ(erow(worst),ecol(worst))), ...
        full(JJfd(erow(worst),ecol(worst))));
    col=ecol(worst); xp=XX; delta=p.solver.fd*max(abs(XX(col)),GG.xs(col));
    xp(col)=xp(col)+delta;
    alone=(residual(xp,XXold,time,dt,GG,BB,p,true)-RR)/delta*GG.xs(col);
    fprintf('Uncompressed FD %.8g; mixed colour contamination %.8g.\n',alone(erow(worst)), ...
        full(JJfd(erow(worst),col))-alone(erow(worst)));
end
assert(freezeJacerr<1e-5,'Frozen-joint tangent differs from independent coloured FD.');
aa=state(TT,hh,GG,p); [~,~,QT]=liningLeak(hh,aa,BB,p);
remoteColumns=find(sum(abs(QT),1)>0);
assert(any(~ismember(remoteColumns,GG.ids)), ...
    'Frozen-joint derivative omitted concrete temperature DOFs.');
report=struct('passed',true,'jacobianRelativeError',jacerr,'sealedWaterDrift_m2',sum(a.W)-W0, ...
    'maxStepWaterDefect_m2',maxW,'maxStepEnergyDefect_J_per_m',maxE,'sealedMaxIce',max(a.thetaI), ...
    'sectorLengthSpread_m',max(segmentLengths)-min(segmentLengths), ...
    'finiteWidthJointFlow_m2_per_s',ledger.byJoint,'acceptedDtRecovered_day',recovered/p.day, ...
    'realMeshJacobianRelativeError',realMeshJacerr, ...
    'frozenJointJacobianRelativeError',freezeJacerr);
% 12 A common enthalpy reference shift must only add the mass residual.
% This detects incompatible control volumes even when global budgets close.
p=defaults(); p.thm.enabled=false; p.hydroFar=false; p.topSeepage=false; p.tunnelDrain=false;
p.geothermalFlux=0; p.mechanics.enable=false;
G=geometry(SmallSoilMesh(),p); B=boundaries(G,SmallSoilMesh(),p);
X=[-.2;.3;-.1;.2;-3;-2;-2;-1]; old=X;
old(1:4)=old(1:4)+.05; old(5:8)=old(5:8)+.1;
[r0,~]=residual(X,old,300,300,G,B,p,false);
p.enthalpyReference_Jkg=1.7e5;
[r1,~]=residual(X,old,300,300,G,B,p,false);
datumDefect=(r1(1:4)-r0(1:4))-p.rho*p.enthalpyReference_Jkg*r0(5:8);
report.enthalpyDatumPatchRelativeError=norm(datumDefect,inf)/max(1,norm(r0(1:4),inf));
assert(report.enthalpyDatumPatchRelativeError<1e-10,'Incompatible heat and water control-volume mapping.');
disp(report);
end

function S=SmallSoilMesh()
S=struct('pAll',[0 0;1 0;0 1;1 1],'tSoil',[1 2 3;2 4 3]);
end

function s=thmDefaults()
% Numerical example; no gasket failure/irreversible freeze-thaw calibration.
s=struct('enabled',true,'phasePressure','bishopIce', ...
    'contactNormalStiffness_Pa_per_m',1e11,'initialJointGap_m',0, ...
    'boltNormalStiffness_N_per_m2',1e9,'boltPreload_N_per_m',2e4, ...
    'boltShearStiffness_N_per_m2',0,'boltRadialFraction',.75, ...
    'openingTransmissionFraction',1e-3,'slitPhaseTransitionWidth_C',.05, ...
    'maxOmittedSlitVolumeFraction',.01, ...
    'forceResidualScale',1e-6,'displacementScale_m',.001, ...
    'maximumVolumetricStrain',.05,'maximumPrincipalStrain',.05, ...
    'maximumDisplacementGradient',.1,'permeabilityPorosityFeedback',false, ...
    'includeSelfWeight',true);
end

function G=thmGeometry(G,S,p)
assert(numel(p.leakage.jointAnglesRad)==8, ...
    'THM:EightJoints','Coupled segment mechanics needs eight longitudinal joints.');
[M,G.mechanicalMeshReport]=build_segment_mechanical_mesh(S,p.leakage.jointAnglesRad);
rs=struct('E_Pa',p.mechanics.Econc,'poissonRatio',p.mechanics.nuConc, ...
    'thermalExpansion_per_C',p.mechanics.alphaConc, ...
    'referenceTemperatureGlobal_C',zeros(G.N,1), ...
    'foundationNormalStiffness_Pa_per_m',p.mechanics.kn, ...
    'foundationTangentialStiffness_Pa_per_m',p.mechanics.kt, ...
    'contactNormalStiffness_Pa_per_m',p.thm.contactNormalStiffness_Pa_per_m, ...
    'initialJointGap_m',p.thm.initialJointGap_m, ...
    'boltNormalStiffness_N_per_m2',p.thm.boltNormalStiffness_N_per_m2, ...
    'boltPreload_N_per_m',p.thm.boltPreload_N_per_m, ...
    'boltShearStiffness_N_per_m2',p.thm.boltShearStiffness_N_per_m2, ...
    'boltRadialFraction',p.thm.boltRadialFraction);
G.ring=solve_segment_contact_fem('assemble',M,rs);
G.gapProfiles=build_joint_gap_profiles(G.ring,G.jointThermalPaths);
[~,~,channelInventory]=joint_channel_inventory( ...
    zeros(G.gapProfiles.nj,G.gapProfiles.nq),G.jointLaw.jointFreeze, ...
    G.jointLaw.jointStorage,G.jointLaw.jointConductance,p.liningThickness, ...
    p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
G.baseChannelWaterVolume=channelInventory.waterVolume_m2;
rhoBulk=(1-p.porosity)*p.rhoS+p.porosity*p.rho;
ss=struct('E',p.mechanics.Esoil,'nu',p.mechanics.nuSoil, ...
    'alphaThermal',p.mechanics.alphaSoil,'biot',p.mechanics.biot, ...
    'rhoBulk',rhoBulk,'gravity',[0;p.g]*p.thm.includeSelfWeight, ...
    'phasePressure',p.thm.phasePressure,'phaseParameters',p);
G.msoil=soil_poroelastic_fem(G.xy(G.ids,:),G.local, ...
    20*ones(G.ns,1),ones(G.ns,1),20*ones(G.ns,1),ss);
[~,G.strainOperator]=soil_small_strain_diagnostics(G.xy(G.ids,:),G.local,zeros(G.ns,2));
G.umap=sparse(udofs(G.ids)',(1:2*G.ns)',1,2*G.N,2*G.ns);
G.nuSoil=2*G.ns; G.nuRing=2*size(M.p,1); G.nm=G.nuSoil+G.nuRing;
G.uIndex=(G.N+G.ns+(1:G.nm))';
G.referenceU=zeros(G.nm,1); G.referenceT=zeros(G.N,1);
ringVolume=zeros(size(M.p,1),1);
for e=1:size(M.q,1)
    ids=M.q(e,:); [~,v]=quadThermal(M.p(ids,:));
    ringVolume(ids)=ringVolume(ids)+v;
end
G.ringBody=zeros(G.nuRing,1);
G.ringBody(2:2:end)=ringVolume*p.rhoC*p.g*p.thm.includeSelfWeight;
xy=G.xy(G.ids,:); tol=1e-8;
sides=find(abs(xy(:,1)-min(xy(:,1)))<tol|abs(xy(:,1)-max(xy(:,1)))<tol);
bottom=find(abs(xy(:,2)-max(xy(:,2)))<tol);
G.mfixed=unique([2*sides-1;2*bottom-1;2*bottom]);
G.mfree=setdiff((1:G.nm)',G.mfixed);
G.mKss=G.msoil.K+G.umap'*G.ring.KsoilSoil*G.umap;
G.mKrs=G.ring.KringSoil*G.umap;
G.qElementSelect=sparse((1:size(G.scv.P,1))',G.scv.element,1,size(G.scv.P,1),numel(G.area));
G.poreDivergence=G.qElementSelect*G.msoil.V;
nr=size(M.jointNodePairs,1); nj=8;
G.gapJointMap=sparse(repelem((1:nj)',nr),(1:nr*nj)',1,nj,nr*nj)*spones(G.ring.gapNodalMap);
G.gapJointMap=spones(G.gapJointMap);
% Include every potentially active contact/bolt, not only the initial branch.
R=G.ring;
ringPattern=spones(R.Kelastic+R.KringRing)+spones(R.gapGaussMap'*R.gapGaussMap) ...
    +spones(R.boltNormalMap'*R.boltNormalMap)+spones(R.boltTangentialMap'*R.boltTangentialMap);
Pm=[spones(G.mKss),spones(G.mKrs');spones(G.mKrs),ringPattern];
soilLoad=spones(spones(G.msoil.A_pressure)+spones(G.msoil.A_temperature))*spones(G.msoil.elementMeanNodal);
LmT=[spones(soilLoad)*G.soilSelect;spones(R.thermalLoadMap)];
LmH=[spones(G.msoil.A_pressure)*spones(G.msoil.elementMeanNodal);sparse(G.nuRing,G.ns)];
Pvol=spones(spones(G.scv.deposit)*spones(G.poreDivergence));
jointSoil=spones(G.jointThermalPaths.interpolation); %#ok<NASGU>
% Every joint path temperature/entry heat row may depend on its opening.
Ppath=G.jointThermalPaths.interpolation; nq=G.jointThermalPaths.nq;
jointThermal=spones(kron(speye(nj),ones(1,nq))*spones(Ppath));
Jsoil=spones(buildJointSoilPattern(G,S,p));
Pring=[spones(G.soilSelect'*Jsoil+jointThermal');Jsoil]*G.gapJointMap;
Pmu=[G.soilSelect'*Pvol;Pvol];
top=[G.pattern,spones([Pmu,Pring])];
bottom=[spones([LmT,LmH]),spones(Pm)];
G.pattern=spones([top;bottom]+speye(G.N+G.ns+G.nm));
G.rs=[G.rs;ones(G.nm,1)];
G.xs=[G.xs;p.thm.displacementScale_m*ones(G.nm,1)];
G.groups=colorColumns(G.pattern); [G.jrow,G.jcol]=find(G.pattern);
G.jnz=cell(size(G.groups)); G.jloc=cell(size(G.groups));
for c=1:numel(G.groups)
    [tf,loc]=ismember(G.jcol,G.groups{c}); G.jnz{c}=find(tf); G.jloc{c}=loc(tf);
end
end

function W=buildJointSoilPattern(G,S,p)
Q=build_joint_drain_quadrature(S,G.ids,p.leakage, ...
    sort(mod(p.leakage.jointAnglesRad(:),2*pi)));
nj=8; nqp=numel(Q.weight);
J=sparse((1:nqp)',Q.jointID,1,nqp,nj);
W=spones(Q.phi'*J);
end

function a=thmState(a,T,h,U,G,p,wantDerivatives)
us=U(1:G.nuSoil); ur=U(G.nuSoil+1:end);
increment=us-G.referenceU(1:G.nuSoil);
epsElement=G.msoil.V*increment; epsQ=G.qElementSelect*epsElement;
Vp=p.porosity+p.mechanics.biot*epsQ; J=1+epsQ; phi=Vp./J;
assert(all(Vp>0&Vp<J)&&max(abs(epsQ))<=p.thm.maximumVolumetricStrain, ...
    'THM:StrainDomain','Small-strain pore-volume domain exceeded.');
a.strainDiagnostics=soil_small_strain_diagnostics(G.strainOperator,increment);
a.absoluteStrainDiagnostics=soil_small_strain_diagnostics(G.strainOperator,us);
assert(a.strainDiagnostics.maxAbsPrincipalStrain<=p.thm.maximumPrincipalStrain ...
    &&a.strainDiagnostics.maxGradientNorm<=p.thm.maximumDisplacementGradient ...
    &&a.absoluteStrainDiagnostics.maxAbsPrincipalStrain<=p.thm.maximumPrincipalStrain ...
    &&a.absoluteStrainDiagnostics.maxGradientNorm<=p.thm.maximumDisplacementGradient, ...
    'THM:StrainDomain','Full strain or displacement-gradient domain exceeded; a finite-deformation formulation is required.');
q=thm_bulk_constitutive(G.scv.P*T(G.ids),G.scv.P*h,epsQ,p,a.quadrature,wantDerivatives);
a.quadrature=q; a.referencePoreVolume=Vp; a.currentPorosityQuadrature=phi;
a.elementVolumetricStrain=epsElement;
epsN=G.msoil.nodeProjection*epsElement;
phiN=(p.porosity+p.mechanics.biot*epsN)./(1+epsN);
a.thetaW=phiN.*a.sl; a.thetaI=phiN.*a.si; a.thetaA=phiN.*a.sg;
a.total=a.thetaW+a.thetaI; a.currentPorosity=phiN;
state=struct('temperatureGlobal_C',T,'referenceTemperatureGlobal_C',G.referenceT, ...
    'ringDisplacement_m',ur,'soilDisplacementGlobal_m',G.umap*us, ...
    'additionalRingLoad_N_per_m',G.ringBody);
a.contact=solve_segment_contact_fem('evaluate',G.ring,state);
paths=G.jointThermalPaths; nq=paths.nq; nj=paths.nj;
profile=reshape(paths.interpolation*T,nq,nj)';
wholePathOpen=all(a.contact.jointGap_m>0,1)';
[factor,dfactor]=joint_freeze_conductance(profile,G.jointLaw.jointFreeze);
options=struct('openingTransmissionFraction',p.thm.openingTransmissionFraction, ...
    'pathLength_m',G.gapProfiles.actualPathLength_m, ...
    'sampleWeights',G.gapProfiles.sampleWeights,'openPathMask',wholePathOpen, ...
    'waterDensity_kg_m3',p.rho,'gravity_m_s2',p.g, ...
    'referenceTemperature_C',p.leakage.jointFreeze.referenceTemperature_C, ...
    'freezeEnabled',p.leakage.jointFreeze.iceEnabled, ...
    'viscosityEnabled',p.leakage.jointFreeze.viscosityEnabled, ...
    'bulkFreezingTemperature_C',p.leakage.jointFreeze.bulkFreezingTemperature_C, ...
    'iceDensity_kg_m3',p.leakage.jointFreeze.iceDensity_kg_m3, ...
    'latentHeat_J_kg',p.leakage.jointFreeze.latentHeat_J_kg, ...
    'meltingTemperature_K',p.leakage.jointFreeze.meltingTemperature_K, ...
    'transitionWidth_C',p.thm.slitPhaseTransitionWidth_C,'baseMobilityFactor',factor, ...
    'baseMobilityDerivative_dT',dfactor);
[a.effectiveJointConductance,a.jointConductanceU,dCdT,a.openingHydraulics]= ...
    joint_opening_p1_hydraulics(G.jointLaw.jointConductance,G.gapProfiles,ur,profile,options);
toJoint=@(v) sparse(repelem((1:nj)',nq),(1:nj*nq)',reshape(v',[],1),nj,nj*nq);
a.jointConductanceT=toJoint(dCdT)*paths.interpolation;
% Additional open slit is a zero-storage boundary approximation. Enforce a
% quantitative volume ceiling instead of silently omitting significant water.
[volume,~,a.openingGeometry]=joint_profile_geometric_volume( ...
    G.gapProfiles,ur,G.jointLaw.jointConductance,p.thm.openingTransmissionFraction);
a.exactOpenSlitGeometricVolume_m2=volume;
baseVolume=G.baseChannelWaterVolume;
ratio=zeros(nj,1); active=baseVolume>0; ratio(active)=volume(active)./baseVolume(active);
a.maxOmittedSlitVolumeFraction=max(ratio);
assert(a.maxOmittedSlitVolumeFraction<=p.thm.maxOmittedSlitVolumeFraction, ...
    'THM:SlitStorageDomain','New open-channel volume is too large for the zero-storage boundary approximation; use a pressure/storage interface model.');
end

function [r,J]=thmMechanicalResidual(T,h,U,G,p,a)
us=U(1:G.nuSoil); ur=U(G.nuSoil+1:end);
if nargin<6
    aq=frozen_soil_constitutive(G.scv.P*T(G.ids),G.scv.P*h,p);
    a=struct('quadrature',aq);
end
aq=a.quadrature;
Tq=G.scv.P*T(G.ids); hq=G.scv.P*h;
pw=p.rho*p.g*hq; pi=pw-p.rho*p.L/p.Tm*min(Tq,0);
if strcmp(p.thm.phasePressure,'bishopIce')
    pressure=aq.sl.*pw+aq.si.*pi;
else
    pressure=pw;
end
fsoil=G.msoil.fbody+G.msoil.A_temperature* ...
    (G.msoil.elementMeanNodal*(T(G.ids)-G.referenceT(G.ids))) ...
    +G.msoil.A_pressure*(G.scv.elementAverage*pressure);
if isfield(a,'contact')
    e=a.contact;
else
    state=struct('temperatureGlobal_C',T,'referenceTemperatureGlobal_C',G.referenceT, ...
        'ringDisplacement_m',ur,'soilDisplacementGlobal_m',G.umap*us, ...
        'additionalRingLoad_N_per_m',G.ringBody);
    e=solve_segment_contact_fem('evaluate',G.ring,state);
end
r=[G.mKss*us+G.mKrs'*ur-fsoil;e.ringResidual_N_per_m];
if nargout>1
    J=[G.mKss,G.mKrs';G.mKrs,e.ringTangent_N_per_m2];
end
end

function U=thmInitialEquilibrium(T,h,G,p)
U=zeros(G.nm,1); free=G.mfree;
fixedPhase=struct('quadrature',frozen_soil_constitutive(G.scv.P*T(G.ids),G.scv.P*h,p,false));
[r,J]=thmMechanicalResidual(T,h,U,G,p,fixedPhase);
scale=max(1,norm(r(free),inf)); converged=false;
for k=1:60
    if norm(r(free),inf)<1e-10*scale+1e-5, converged=true; break; end
    du=zeros(G.nm,1); du(free)=-J(free,free)\r(free);
    accepted=false; alpha=1;
    for ls=1:24
        candidate=U+alpha*du;
        [rc,Jc]=thmMechanicalResidual(T,h,candidate,G,p,fixedPhase);
        if norm(rc(free),2)<(1-1e-4*alpha)*norm(r(free),2)
            accepted=true; break
        end
        alpha=alpha/2;
    end
    assert(accepted,'THM:InitialEquilibrium','Installed contact equilibrium line search failed.');
    U=candidate; r=rc; J=Jc;
end
assert(converged,'THM:InitialEquilibrium','Installed selfweight/hydraulic/contact equilibrium did not converge.');
end

function ok=thmAdmissible(X,G,p)
ok=all(isfinite(X)); if ~ok||~p.thm.enabled, return; end
U=X(G.uIndex); eps=G.msoil.V*(U(1:G.nuSoil)-G.referenceU(1:G.nuSoil));
Vp=p.porosity+p.mechanics.biot*eps;
ok=all(abs(eps)<=p.thm.maximumVolumetricStrain&Vp>0&Vp<1+eps);
if ~ok, return; end
d=soil_small_strain_diagnostics(G.strainOperator,U(1:G.nuSoil)-G.referenceU(1:G.nuSoil));
ok=d.maxAbsPrincipalStrain<=p.thm.maximumPrincipalStrain ...
    &&d.maxGradientNorm<=p.thm.maximumDisplacementGradient;
if ~ok, return; end
sa=soil_small_strain_diagnostics(G.strainOperator,U(1:G.nuSoil));
ok=sa.maxAbsPrincipalStrain<=p.thm.maximumPrincipalStrain ...
    &&sa.maxGradientNorm<=p.thm.maximumDisplacementGradient;
if ~ok, return; end
ur=U(G.nuSoil+1:end); B=G.gapProfiles;
volume=joint_profile_geometric_volume(B,ur, ...
    G.jointLaw.jointConductance,p.thm.openingTransmissionFraction);
active=G.baseChannelWaterVolume>0;
ok=all(volume(active)<=p.thm.maxOmittedSlitVolumeFraction*G.baseChannelWaterVolume(active));
end

function [xp,delta]=thmPerturbation(X,cols,G,p)
% Normally all positive perturbations are feasible. Close to a constitutive
% boundary choose a feasible one-sided direction per column, then reduce the
% coloured increment if their combination exceeds a joint-volume/strain
% bound. An accepted state never bypasses the physical-domain checks.
delta=p.solver.fd*max(abs(X(cols)),G.xs(cols));
xp=X; xp(cols)=xp(cols)+delta;
if thmAdmissible(xp,G,p), return; end
for k=1:numel(cols)
    d=delta(k); feasible=false;
    for attempt=1:24
        probe=X; probe(cols(k))=probe(cols(k))+d;
        if thmAdmissible(probe,G,p), feasible=true; break; end
        probe(cols(k))=X(cols(k))-d;
        if thmAdmissible(probe,G,p), d=-d; feasible=true; break; end
        d=d/2;
    end
    assert(feasible,'THM:TangentDomain','No admissible one-sided perturbation at the current constitutive boundary.');
    delta(k)=d;
end
for attempt=1:24
    xp=X; xp(cols)=xp(cols)+delta;
    if thmAdmissible(xp,G,p), return; end
    delta=delta/2;
end
error('THM:TangentDomain','Coloured constitutive perturbations exceed the admissible domain.');
end

function m=thmMechanicalOutput(X,G,p)
T=X(1:G.N); h=X(G.N+(1:G.ns)); U=X(G.uIndex);
a=state(T,h,G,p,U); r=thmMechanicalResidual(T,h,U,G,p,a);
stressState=struct('temperatureGlobal_C',T,'referenceTemperatureGlobal_C',G.referenceT, ...
    'ringDisplacement_m',U(G.nuSoil+1:end),'soilDisplacementGlobal_m',G.umap*U(1:G.nuSoil), ...
    'additionalRingLoad_N_per_m',G.ringBody,'computeStress',true);
a.contact=solve_segment_contact_fem('evaluate',G.ring,stressState);
m=struct('formulation','Monolithic equal-density small-strain THM with eight contact blocks', ...
    'soilU_m',reshape(U(1:G.nuSoil),2,[])', ...
    'soilIncrementU_m',reshape(U(1:G.nuSoil)-G.referenceU(1:G.nuSoil),2,[])', ...
    'ringU_m',reshape(U(G.nuSoil+1:end),2,[])', ...
    'referenceU',G.referenceU,'referenceTemperature',G.referenceT, ...
    'soilElementVolumetricStrain',a.elementVolumetricStrain, ...
    'soilStrainDiagnostics',a.strainDiagnostics, ...
    'absoluteSoilStrainDiagnostics',a.absoluteStrainDiagnostics, ...
    'installedSoilStrainDiagnostics',G.referenceStrainDiagnostics, ...
    'currentSoilPorosity',a.currentPorosity, ...
    'contact',a.contact,'openingHydraulics',a.openingHydraulics, ...
    'openingGeometry',a.openingGeometry, ...
    'exactOpenSlitGeometricVolume_m2',a.exactOpenSlitGeometricVolume_m2, ...
    'mesh',G.ring.mesh,'maximumFreeResidual_N_per_m',max(abs(r(G.mfree))), ...
    'omittedNewChannelVolumeFraction',a.maxOmittedSlitVolumeFraction, ...
    'status','Example contact/bolt parameters; no gasket failure calibration or construction sequence', ...
    'limitations','No segregated ice lenses, plasticity, fracture, density-expansion model, mechanical heating or resolved open-channel storage.');
end

function report=thmSelftest()
S=make_thm_verification_mesh();
o=joint_supply_th_model(struct('mesh',S,'durationDays',.05,'dtInitial',.01*86400, ...
    'dtMax',.01*86400,'outputEveryDays',.05,'plot',false,'saveOutput',false,'verbose',false, ...
    'waterTableDepth',-2,'topSeepage',false,'geothermalFlux',0, ...
    'topMin',.2,'topMax',.2,'tunnelMin',-.1,'tunnelMax',-.1, ...
    'initial',struct('surfaceT',-.03,'gradient',.01),'mechanics',struct('enable',false)));
out=TH_2D_segment_v16(o); p=out.parameters;
off=o; off.thm=struct('enabled',false); uncoupled=TH_2D_segment_v16(off);
G=geometry(S,p); B=boundaries(G,S,p); G=jointCouplingPattern(G,B);
G.jointThermalPaths=B.jointThermalPaths;
Q=B.jointQuadrature; nj=8; nq=B.jointThermalPaths.nq;
Jmap=sparse((1:numel(Q.weight))',Q.jointID,1,numel(Q.weight),nj);
[~,~,inventory]=joint_channel_inventory(zeros(nj,nq),p.leakage.jointFreeze, ...
    p.leakage.jointStorage,p.leakage.jointConductance,p.liningThickness,p.rho,p.L,p.cw,p.ci,p.rhoC*p.cc,p.g);
G.channelWaterNodeMass=Q.phi'*spdiags(Q.weight,0,numel(Q.weight),numel(Q.weight))*Jmap*inventory.waterVolume_m2;
G=thmGeometry(G,S,p);
G.referenceU=out.mechanics.referenceU; G.referenceT=out.mechanics.referenceTemperature;
U=[reshape(out.mechanics.soilU_m',[],1);reshape(out.mechanics.ringU_m',[],1)];
X=[out.T;out.h;U]; old=X;
X(G.N+(1:G.ns))=X(G.N+(1:G.ns))+.05*sin(G.xy(G.ids,1));
% Activate a genuinely open defect, so the independent tangent check also
% covers ring U -> aperture -> flow and remote heat-row dependencies.
X(1:G.N)=2;
pair=G.ring.mesh.jointNodePairs(:,:,7);
normal=[-sin(p.leakage.jointAnglesRad(7));cos(p.leakage.jointAnglesRad(7))];
trialU=reshape(X(G.N+G.ns+G.nuSoil+(1:G.nuRing)),2,[])';
for k=1:size(pair,1)
    meanU=(trialU(pair(k,1),:)+trialU(pair(k,2),:))/2;
    trialU(pair(k,1),:)=meanU-.5*6e-4*normal';
    trialU(pair(k,2),:)=meanU+.5*6e-4*normal';
end
X(G.N+G.ns+G.nuSoil+(1:G.nuRing))=reshape(trialU',[],1);
old(1:G.N)=old(1:G.N)+.001;
dt=864; time=4320;
[R,D]=residual(X,old,time,dt,G,B,p,true);
trialState=state(X(1:G.N),X(G.N+(1:G.ns)),G,p,X(G.uIndex));
assert(trialState.openingHydraulics.openConductance(7)>0);
assert(all(trialState.effectiveJointConductance([1:6,8])==0));
coldX=X; coldX(1:G.N)=-1;
pNoIce=p; pNoIce.leakage.jointFreeze.iceEnabled=false;
pNoIce.leakage.jointFreeze.viscosityEnabled=false;
coldState=state(coldX(1:G.N),coldX(G.N+(1:G.ns)),G,pNoIce,coldX(G.uIndex));
iceSwitchError=abs(coldState.openingHydraulics.openConductance(7) ...
    /coldState.openingHydraulics.referenceOpenConductance(7)-1);
assert(iceSwitchError<1e-13,'THM:FreezeSwitch','Disabled ice/viscosity still suppresses the open slit.');
% Both defect and new-slit branches must use the same prescribed bulk
% freezing temperature. Below zero does not necessarily mean frozen.
pDepressed=p;
pDepressed.leakage.jointFreeze.bulkFreezingTemperature_C=-.4;
pDepressed.leakage.jointFreeze.viscosityEnabled=false;
mildState=state(-.2*ones(G.N,1),X(G.N+(1:G.ns)),G,pDepressed,X(G.uIndex));
subzeroSlitMobility=mildState.openingHydraulics.openConductance(7) ...
    /mildState.openingHydraulics.referenceOpenConductance(7);
assert(abs(subzeroSlitMobility-1)<1e-13,'THM:JointBulkFreezing', ...
    'A slit above the prescribed depressed freezing temperature was blocked.');
deepColdState=state(-.6*ones(G.N,1),X(G.N+(1:G.ns)),G,pDepressed,X(G.uIndex));
assert(deepColdState.openingHydraulics.openConductance(7)==0,'THM:JointBulkFreezing', ...
    'A slit below the prescribed freezing temperature still conducts.');
% An overlarge trial slit must be backtracked before evaluating a residual.
overlarge=X; largeU=trialU;
for k=1:size(pair,1)
    meanU=(largeU(pair(k,1),:)+largeU(pair(k,2),:))/2;
    largeU(pair(k,1),:)=meanU-.5*4e-3*normal';
    largeU(pair(k,2),:)=meanU+.5*4e-3*normal';
end
overlarge(G.N+G.ns+G.nuSoil+(1:G.nuRing))=reshape(largeU',[],1);
assert(~thmAdmissible(overlarge,G,p),'THM:TrialDomain','Oversized trial aperture was not rejected.');
J=jacobian(X,old,time,dt,G,B,p,R);
% Entire small-system uncompressed FD checks both colouring and sparsity,
% including remote ring openings and adjacent soil-strain storage columns.
Jfull=zeros(numel(X),numel(X));
for j=1:numel(X)
    [xp,delta]=thmPerturbation(X,j,G,p);
    Jfull(:,j)=(residual(xp,old,time,dt,G,B,p,true)-R)/delta*G.xs(j);
end
jacError=norm(full(J)-Jfull,inf)/max(1,norm(Jfull,inf));
assert(jacError<1e-5,'THM:Jacobian','Monolithic tangent/pattern misses coupled dependencies.');
% Water/ice reference energy must remain conjugate to mass when pores move.
r0=residual(X,old,time,dt,G,B,p,false);
p.enthalpyReference_Jkg=1.7e5; r1=residual(X,old,time,dt,G,B,p,false);
freeSoil=setdiff((1:G.ns)',B.far);
patch=(r1(G.ids(freeSoil))-r0(G.ids(freeSoil))) ...
    -p.rho*p.enthalpyReference_Jkg*r0(G.N+freeSoil);
datumError=norm(patch,inf)/max(1,norm(r0(1:G.N),inf));
assert(datumError<1e-9,'THM:Datum','Mechanical pore storage breaks energy-reference invariance.');
headFeedback=max(abs(out.h-uncoupled.h));
assert(headFeedback>1e-7,'THM:Feedback','Mechanical deformation has no effect on water solution.');
assert(max(abs(out.state.currentPorosity-p.porosity))>1e-9);
assert(max(abs(out.history(:,8)))<1e-7&&max(abs(out.history(:,9)))<10);
assert(out.mechanics.maximumFreeResidual_N_per_m<.01);
report=struct('passed',true,'jacobianRelativeError',jacError,'colourCount',numel(G.groups), ...
    'enthalpyDatumRelativeError',datumError,'maximumHeadFeedback_m',headFeedback, ...
    'maxPorosityChange',max(abs(out.state.currentPorosity-p.porosity)), ...
    'maximumWaterDefect_m2',max(abs(out.history(:,8))), ...
    'maximumEnergyDefect_J_per_m',max(abs(out.history(:,9))), ...
    'maximumFreeMechanicalResidual_Nm',out.mechanics.maximumFreeResidual_N_per_m, ...
    'maxVolumetricStrain',max(abs(out.state.elementVolumetricStrain)), ...
    'trialOpenConductance_m_s',trialState.openingHydraulics.openConductance(7), ...
    'oversizedTrialRejected',true, ...
    'disabledIceViscositySwitchError',iceSwitchError, ...
    'depressedFreezingSubzeroSlitMobility',subzeroSlitMobility, ...
    'depressedFreezingColdSlitConductance_m_s',deepColdState.openingHydraulics.openConductance(7), ...
    'coupledOutput',out,'uncoupledOutput',uncoupled, ...
    'localWaterDefectAtTrial_m2',max(abs(D.localWaterDefect)));
disp(rmfield(report,{'coupledOutput','uncoupledOutput'}));
end

