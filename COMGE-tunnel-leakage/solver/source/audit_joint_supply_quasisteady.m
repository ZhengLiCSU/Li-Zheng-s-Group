function audit=audit_joint_supply_quasisteady(studyFolder,caseNames)
%AUDIT_JOINT_SUPPLY_QUASISTEADY Conditional fixed-mobility water audit.
% Read complete v14 twenty-year cases and their twelve monthly snapshots.
% At each snapshot hold T, nodal soil K(T,hTH), joint thermal factor and far
% water table fixed. Solve a no-storage steady head equation with the same
% triangle-mean K operator and consistent unilateral joint/top drainage.
% This is not an independent Richards/TH calculation: K and temperatures
% come from TH and are NOT updated as the conditional steady head changes.
% Differences therefore characterize removal of transient storage at fixed
% instantaneous mobility. Twelve snapshots cannot give an exact annual
% quasi-steady leakage integral; no annual integral is reported.
%
% audit=audit_joint_supply_quasisteady(STUDYFOLDER,{'base','supplySeasonal'});
% Another fine-mesh case may be named if it has its own full v14 history.
root=fileparts(mfilename('fullpath')); addpath(root);
if nargin<1||isempty(studyFolder)
    studyFolder=fullfile(root,'joint_supply_study_20260930');
end
if nargin<2||isempty(caseNames), caseNames={'base','supplySeasonal'}; end
studyFolder=char(studyFolder); caseNames=cellstr(string(caseNames));
assert(~isempty(caseNames)&&numel(unique(caseNames))==numel(caseNames), ...
    'Supply distinct complete-history case names.');
outputFolder=fullfile(studyFolder,'fixed_mobility_quasisteady_audit');
if ~isfolder(outputFolder), mkdir(outputFolder); end
rows=struct([]); details=struct([]); record=0;
for caseIndex=1:numel(caseNames)
    name=caseNames{caseIndex};
    sourceFile=fullfile(studyFolder,name,'TH_segment_spinup_results.mat');
    assert(isfile(sourceFile),'Missing complete source: %s.',sourceFile);
    loaded=load(sourceFile,'result'); r=loaded.result;
    assert(r.completedYears==20&&strcmp(r.lastYear.version,'TH_v14_20260930_joint_supply'), ...
        'Each audit source must complete its own twenty-year v14 history.');
    out=r.lastYear; p=out.parameters; S=out.mesh;
    assert(p.hydroFar&&p.tunnelDrain&&p.liningK==0, ...
        'This audit requires prescribed far heads and joint-only tunnel extraction.');
    assert(numel(out.snapshots)==12,'The complete twelve-snapshot twentieth year is required.');
    [B,nodes,far,free,local]=prepareBoundary(S,p);
    assert(isequal(nodes,out.soilNodes(:)), ...
        'Snapshot hydraulic node order differs from its own source mesh.');
    phaseLevel=NaN;
    if isfield(S,'jointRefinementReport'), phaseLevel=S.jointRefinementReport.level; end
    for snapshot=1:12
        snap=out.snapshots(snapshot); T=snap.T(:); h=snap.h(:);
        day=snap.timeDay; phaseDay=day-19*365;
        assert(day>19*365&&day<=20*365+1e-7&&numel(h)==numel(nodes), ...
            'Invalid twentieth-year snapshot time or head dimensions.');
        K=frozenSoilMobility(T(nodes),h,p);
        [A,trans]=soilOperator(S.pAll(nodes,:),local,K);
        depth=waterDepth(day,p);
        href=B.soilDepth-depth;
        drawdownActual=href-h; drawdownActual(far)=0;
        [qActual,actualLedger]=evaluate_joint_drain(h,T,B,p,depth);
        topCoeff=B.topW*p.topLeakConductance*double(p.topSeepage);
        topActual=topCoeff.*max(h,0);
        actualNoStorageResidual=-A*drawdownActual+qActual+topActual;
        nearest=find(abs(out.history(:,1)-day)<1e-7,1);
        assert(~isempty(nearest),'Snapshot has no matching accepted TH step.');
        historyLeak=out.history(nearest,10);
        assert(abs(historyLeak-sum(qActual))<=1e-12+1e-8*abs(historyLeak), ...
            'Reconstructed snapshot leakage differs from its source TH ledger.');
        [steadyHead,solver]=steadyFrozenMobility(A,href,drawdownActual,far,free, ...
            T,B,p,depth,topCoeff);
        qSteady=nan(size(qActual)); steadyLedger=[]; massDefect=NaN;
        steadyRate=NaN; headSteady=NaN; farInflow=NaN;
        if solver.converged
            [qSteady,steadyLedger]=evaluate_joint_drain(steadyHead,T,B,p,depth);
            steadyDrawdown=href-steadyHead;
            topSteady=topCoeff.*max(steadyHead,0);
            waterResidual=-A*steadyDrawdown+qSteady+topSteady;
            reaction=zeros(numel(h),1); reaction(far)=-waterResidual(far);
            massDefect=sum(qSteady+topSteady+reaction);
            steadyRate=sum(qSteady); farInflow=-sum(reaction);
            headSteady=steadyLedger.drivingHead(7);
        end
        actualRate=sum(qActual);
        relativeDifference=NaN;
        if solver.converged
            relativeDifference=100*abs(actualRate-steadyRate) ...
                /max([abs(actualRate),abs(steadyRate),realmin]);
        end
        record=record+1;
        row=struct('caseName',string(name),'meshRefinementLevel',phaseLevel, ...
            'snapshotIndex',snapshot,'absoluteDay',day,'phaseDay',phaseDay, ...
            'waterTableDepth_m',depth,'KsatReference_m_per_s',p.Ksat, ...
            'nodalK_min_m_per_s',min(K),'nodalK_max_m_per_s',max(K), ...
            'minimumTransmissivity_m_per_s',min(trans),'negativeTransmissivityEdges',sum(trans<0), ...
            'joint7_T_factor',actualLedger.conductance(7)/p.leakage.jointConductance(7), ...
            'TH_actualRate_m2_per_s',actualRate, ...
            'fixedMobilityQuasiSteadyRate_m2_per_s',steadyRate, ...
            'TH_actual_L_per_day_m',actualRate*86400*1000, ...
            'fixedMobilityQuasiSteady_L_per_day_m',steadyRate*86400*1000, ...
            'snapshotDifference_percentOfLargerRate',relativeDifference, ...
            'TH_joint7_drivingHead_m',actualLedger.drivingHead(7), ...
            'fixedMobilityQuasiSteady_joint7_drivingHead_m',headSteady, ...
            'noStorageResidualAtTH_maxFree_m2_per_s',norm(actualNoStorageResidual(free),inf), ...
            'quasiSteadyFarInflow_m2_per_s',farInflow,'quasiSteadyMassDefect_m2_per_s',massDefect, ...
            'quasiSteadyConverged',solver.converged,'newtonIterations',solver.iterations, ...
            'quasiSteadyMaxFreeResidual_m2_per_s',solver.maximumResidual, ...
            'quasiSteadyScaledHeadResidual_m',solver.scaledHeadResidual, ...
            'solverReason',string(solver.reason));
        detail=struct('caseName',name,'snapshotIndex',snapshot,'sourceFile',sourceFile, ...
            'soilNodes',nodes,'snapshotTemperature_C',T,'TH_head_m',h, ...
            'frozenMobility_m_per_s',K,'steadyHead_m',steadyHead, ...
            'TH_nodalDrain_m2_per_s',qActual,'steadyNodalDrain_m2_per_s',qSteady, ...
            'TH_jointLedger',actualLedger,'steadyJointLedger',steadyLedger,'solver',solver);
        if record==1, rows=row; details=detail;
        else, rows(record)=row; details(record)=detail; end %#ok<AGROW>
        fprintf('Frozen-mobility audit %s snapshot %d: QTH=%.6g, Qsteady=%.6g L/day/m, %s.\n', ...
            name,snapshot,row.TH_actual_L_per_day_m,row.fixedMobilityQuasiSteady_L_per_day_m,solver.reason);
    end
end
comparison=struct2table(rows);
definition=struct('kind','conditional fixed instantaneous mobility quasi-steady audit', ...
    'sourceSolver','TH_v14_20260930_joint_supply','caseNames',{caseNames}, ...
    'mobilityPolicy',['Soil K is reconstructed once from actual snapshot T/h, ' ...
        'Painter-Karra beta=omega=1, VG-Mualem, ice impedance and optional IAPWS viscosity. ' ...
        'It is held fixed during the no-storage head solve.'], ...
    'drainPolicy','same consistent positive-part joint quadrature and nodal top seepage', ...
    'boundaryPolicy','same snapshot water table, side/bottom head, surface no imposed recharge', ...
    'numericalPolicy','exact source triangle-mean K edge operator; no negative weight or K clipping', ...
    'interpretation',['The fixed-mobility result is conditional on TH-derived temperatures ' ...
        'and conductivities. It is not an independent nonlinear Richards or TH solution. ' ...
        'Snapshot differences support a conditional storage diagnostic only. Twelve ' ...
        'snapshots cannot establish an exact annual leakage or its integration error. ' ...
        'Failed solves remain marked with NaN and are not replaced by hydrostatic heads.']);
audit=struct('definition',definition,'comparison',comparison,'details',details);
writetable(comparison,fullfile(outputFolder,'fixed_mobility_quasisteady_snapshots.csv'));
save(fullfile(outputFolder,'fixed_mobility_quasisteady_audit.mat'),'audit','-v7.3');
disp(comparison(:,{'caseName','snapshotIndex','TH_actual_L_per_day_m', ...
    'fixedMobilityQuasiSteady_L_per_day_m','snapshotDifference_percentOfLargerRate', ...
    'quasiSteadyConverged'}));
end

function [B,nodes,far,free,local]=prepareBoundary(S,p)
xy=S.pAll; nodes=unique(S.tSoil(:)); ns=numel(nodes); N=size(xy,1);
map=zeros(N,1); map(nodes)=1:ns; local=map(S.tSoil);
soilXY=xy(nodes,:); low=min(soilXY,[],1); high=max(soilXY,[],1); tol=1e-8;
far=find(abs(soilXY(:,1)-low(1))<tol|abs(soilXY(:,1)-high(1))<tol ...
    |abs(soilXY(:,2)-high(2))<tol);
free=setdiff((1:ns)',far);
top=find(abs(soilXY(:,2)-low(2))<tol); [~,order]=sort(soilXY(top,1)); top=top(order);
topW=zeros(ns,1); lengths=hypot(diff(soilXY(top,1)),diff(soilXY(top,2)));
topW(top(1:end-1))=topW(top(1:end-1))+lengths/2;
topW(top(2:end))=topW(top(2:end))+lengths/2;
interface=unique(S.soilTunnelCoupleTable(:,2:3)); inner=S.tunnelInnerRingGlobal(:);
outerR=hypot(xy(interface,1)-S.xc,xy(interface,2)-S.yc);
innerR=mean(hypot(xy(inner,1)-S.xc,xy(inner,2)-S.yc));
drop=zeros(ns,1); drop(map(interface))=-(outerR-innerR).*(xy(interface,2)-S.yc)./outerR;
sector=sort(mod(p.leakage.jointAnglesRad(:)-p.leakage.segmentStartAngleRad,2*pi));
sector=mod(sector+p.leakage.segmentStartAngleRad,2*pi);
B=struct('soilSelect',sparse((1:ns)',nodes,1,ns,N), ...
    'soilDepth',soilXY(:,2)-min(xy(:,2)),'liningZDrop',drop,'topW',topW, ...
    'jointQuadrature',build_joint_drain_quadrature(S,nodes,p.leakage,sector));
if p.leakage.jointFreeze.enabled
    B.jointThermalPaths=build_joint_thermal_paths(S,p.leakage.jointAnglesRad, ...
        p.leakage.jointFreeze.quadraturePoints);
end
end

function K=frozenSoilMobility(T,h,p)
m=1-1/p.vgN; sr=p.thetaR/p.porosity;
r=sr+(1-sr)*(1+(p.vgAlpha*max(-h,0)).^p.vgN).^(-m);
hFreeze=-(p.L/(p.g*p.Tm))*max(-T,0);
F=sr+(1-sr)*(1+(p.vgAlpha*max(-hFreeze,0)).^p.vgN).^(-m);
sl=min(r,F); si=max(0,1-F./r);
Se=max(0,min(1,(sl-sr)/(1-sr)));
kr=sqrt(Se).*(1-(1-Se.^(1/m)).^m).^2;
K=p.Ksat*kr.*10.^(-p.iceImpedance*si);
if p.soilViscosity.enabled
    K=K.*(water_viscosity_iapws(p.soilViscosity.referenceTemperature_C) ...
        ./water_viscosity_iapws(T));
end
assert(all(isfinite(K))&&all(K>=0),'Snapshot mobility contains invalid values.');
end

function [A,trans]=soilOperator(xy,tri,K)
a=xy(tri(:,1),:); b=xy(tri(:,2),:); c=xy(tri(:,3),:);
twice=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(c(:,1)-a(:,1)).*(b(:,2)-a(:,2));
assert(all(abs(twice)>1e-12),'Degenerate soil triangle.');
area=abs(twice)/2;
gradient=cat(3,[b(:,2)-c(:,2),c(:,1)-b(:,1)]./twice, ...
    [c(:,2)-a(:,2),a(:,1)-c(:,1)]./twice, ...
    [a(:,2)-b(:,2),b(:,1)-a(:,1)]./twice);
ke=mean(K(tri),2);
pairs=[tri(:,[1 2]);tri(:,[2 3]);tri(:,[3 1])];
value=-[ke.*area.*sum(gradient(:,:,1).*gradient(:,:,2),2); ...
    ke.*area.*sum(gradient(:,:,2).*gradient(:,:,3),2); ...
    ke.*area.*sum(gradient(:,:,3).*gradient(:,:,1),2)];
[edges,~,edgeMap]=unique(sort(pairs,2),'rows');
trans=accumarray(edgeMap,value,[size(edges,1) 1]); nf=size(edges,1); ns=size(xy,1);
inc=sparse([edges(:,1);edges(:,2)],[(1:nf)';(1:nf)'], ...
    [ones(nf,1);-ones(nf,1)],ns,nf);
A=inc*spdiags(trans,0,nf,nf)*inc';
end

function depth=waterDepth(day,p)
depth=p.waterTableDepth;
if ~isempty(p.forcing.waterTableDepth_m)
    phase=mod(day,p.yearDays);
    depth=interp1(p.forcing.day,p.forcing.waterTableDepth_m,phase,'linear');
end
end

function [head,info]=steadyFrozenMobility(A,href,start,far,free,T,B,p,depth,topCoeff)
info=struct('converged',false,'iterations',0,'maximumResidual',NaN, ...
    'scaledHeadResidual',NaN,'reason','iteration limit');
s=start; s(far)=0; head=nan(size(href));
try
    for iteration=1:50
        currentHead=href-s;
        [q,~,DH]=evaluate_joint_drain(currentHead,T,B,p,depth);
        topFlow=topCoeff.*max(currentHead,0);
        residual=A*s-q-topFlow;
        J=A+DH+spdiags(topCoeff.*(currentHead>0),0,numel(s),numel(s));
        diagonal=full(diag(J(free,free)));
        assert(all(isfinite(diagonal))&&all(diagonal>0), ...
            'Frozen-mobility steady matrix has zero/invalid free diagonals.');
        scaleResidual=max(abs(residual(free))./diagonal);
        info.iterations=iteration; info.maximumResidual=norm(residual(free),inf);
        info.scaledHeadResidual=scaleResidual;
        if scaleResidual<1e-7
            info.converged=true; info.reason='converged'; head=currentHead; return;
        end
        d=1./sqrt(diagonal); D=spdiags(d,0,numel(free),numel(free));
        Jscaled=D*J(free,free)*D;
        [R,flag]=chol(Jscaled);
        assert(flag==0,'Frozen-mobility steady matrix failed scaled Cholesky.');
        delta=zeros(size(s)); delta(free)=-d.*(R\(R'\(d.*residual(free))));
        assert(all(isfinite(delta)),'Nonfinite fixed-mobility head increment.');
        objective=energyObjective(s,A,href,T,B,p,depth,topCoeff);
        directional=residual(free)'*delta(free);
        accepted=false;
        for backtrack=0:30
            alpha=2^-backtrack; trial=s+alpha*delta; trial(far)=0;
            merit=energyObjective(trial,A,href,T,B,p,depth,topCoeff);
            if merit<=objective+1e-4*alpha*directional ...
                    +20*eps(max(abs(objective),realmin))
                s=trial; accepted=true; break;
            end
        end
        if ~accepted, info.reason='convex active-set line search failed'; return; end
    end
catch ME
    info.reason=['solve failed: ' ME.message];
end
end

function value=energyObjective(s,A,href,T,B,p,depth,topCoeff)
head=href-s; [~,ledger]=evaluate_joint_drain(head,T,B,p,depth);
Q=B.jointQuadrature;
excess=max(ledger.quadratureHead-ledger.entryHead(Q.jointID),0);
mobility=ledger.conductance(Q.jointID).*Q.weight;
value=.5*(s'*(A*s))+.5*sum(mobility.*excess.^2) ...
    +.5*sum(topCoeff.*max(head,0).^2);
end
