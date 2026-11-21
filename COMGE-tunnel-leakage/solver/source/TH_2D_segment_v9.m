function out = TH_2D_segment_v9(opt)
% TH_2D_SEGMENT_V9 Frozen-ground TH model with total-head lining leakage.
% Usage: out=TH_2D_segment_v9(opt). No PDE/Optimization Toolbox.
% T [degC], h [m LIQUID pressure head], p=rho*g*h [Pa], Z positive down.
% Painter & Karra (2014), Eqs 10-16, beta=omega=1 special case; VG-Mualem.
% Fixed pore volume, equal phase density approximation in the TH equations.
% Mechanics is OPTIONAL ONE-WAY small-strain diagnostic, NOT full THM.
% No independent moisture PDE, pressure clipping, or artificial relaxation.
% Joint conductance and hydraulic head-loss threshold require calibration.
if nargin<1, opt=struct(); end
if ischar(opt) && strcmp(opt,'selftest'), out=selftest(); return; end
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
[X,t0,initInfo]=initialFields(G,B,p);
if p.verbose
    fprintf('TH v9: %d thermal + %d liquid-head DOFs; %d Jacobian colors.\n',G.N,G.ns,numel(G.groups));
    fprintf('Equal-density TH; mechanics=%d (one-way). No cyclic convergence assumed.\n',p.mechanics.enable);
end
start=X; hist=[]; cycles=[]; snaps=struct('timeDay',{},'T',{},'h',{},'thetaW',{},'thetaI',{});
leakHistSegment=zeros(0,p.leakage.segmentCount);
leakHistJoint=zeros(0,numel(p.leakage.jointAnglesRad));
leakHistBody=zeros(0,1);
yearRef=X; cycleConverged=false; t=t0; dt=p.dtInitial; naccepted=0; nrejected=0;
endTime=t0+p.durationDays*p.day; nextReport=t0+p.outputEveryDays*p.day;
nextYear=t0+p.yearDays*p.day; productionEnd=inf;
if ~strcmp(p.initial.mode,'spinup'), productionEnd=endTime; end
while t < min(endTime,productionEnd)-1e-7
    proposedDt=dt;
    target=min([t+dt,endTime,productionEnd,nextReport,nextYear]); d=target-t;
    [Xtrial,info]=step(X,t,d,G,B,p);
    if ~info.ok
        nrejected=nrejected+1; dt=d/2;
        if dt<p.dtMin
            if p.saveOutput
                checkpoint=makeCheckpoint(X,t,G,p,initInfo); %#ok<NASGU>
                save([p.outputPrefix '_FAILED_checkpoint.mat'],'checkpoint','hist','cycles','snaps', ...
                    'leakHistSegment','leakHistJoint','leakHistBody','info');
            end
            error('TH:nonconvergence','Step at day %.5f failed: %s. Last accepted state retained; reduce forcing/inspect parameters.',t/p.day,info.reason);
        end
        if p.verbose, fprintf('  retry day %.4f: dt -> %.5g day (%s)\n',t/p.day,dt/p.day,info.reason); end
        continue
    end
    [~,D]=residual(Xtrial,X,target,d,G,B,p,false);
    X=Xtrial; t=target; naccepted=naccepted+1;
    rec=[t/p.day,d/p.day,info.iter,D.water,D.energy,D.Qout,D.Eout,D.waterError,D.energyError,D.leak,info.norm, ...
        info.jacobianBuilds,info.linearSolves,info.residualEvaluations];
    hist=[hist;rec]; %#ok<AGROW>
    leakHistSegment=[leakHistSegment;D.leakBySegment(:)']; %#ok<AGROW>
    leakHistJoint=[leakHistJoint;D.leakByJoint(:)']; %#ok<AGROW>
    leakHistBody=[leakHistBody;D.leakBody]; %#ok<AGROW>
    % Work-based controller. Use the requested step, not a small remainder
    % imposed by a monthly/yearly output boundary. In v6, an accepted step
    % with moderate effort could remain trapped near dtMin indefinitely.
    effort=max(info.residualEvaluations,1);
    factor=p.solver.dtGrowthSafety*sqrt(p.solver.dtTargetResidualEvals/effort);
    factor=min(p.solver.dtGrowthMax,max(p.solver.dtShrinkMin,factor));
    if info.iter>=ceil(.75*p.solver.maxIter), factor=min(factor,.7); end
    dt=min(p.dtMax,max(p.dtMin,proposedDt*factor));
    if t>=nextReport-1e-7
        a=state(X(1:G.N),X(G.N+1:end),G,p);
        snaps(end+1)=struct('timeDay',t/p.day,'T',X(1:G.N),'h',X(G.N+1:end), ...
            'thetaW',a.thetaW,'thetaI',a.thetaI); %#ok<AGROW>
        if p.verbose
            fprintf('day %8.2f: T=[%7.3f,%7.3f], h=[%8.2f,%8.2f], iceMax=%.4f, leak=%.4g m3/day/m, residual=%.2g\n', ...
                t/p.day,min(X(1:G.N)),max(X(1:G.N)),min(X(G.N+1:end)),max(X(G.N+1:end)),max(a.thetaI),D.leak*p.day,info.norm);
        end
        if p.saveOutput
            checkpoint=makeCheckpoint(X,t,G,p,initInfo); %#ok<NASGU>
            save([p.outputPrefix '_checkpoint.mat'],'checkpoint','hist','cycles','snaps', ...
                'leakHistSegment','leakHistJoint','leakHistBody');
        end
        nextReport=nextReport+p.outputEveryDays*p.day;
    end
    if t>=nextYear-1e-7
        a=state(X(1:G.N),X(G.N+1:end),G,p);
        ar=state(yearRef(1:G.N),yearRef(G.N+1:end),G,p);
        eT=max(abs(X(1:G.N)-yearRef(1:G.N)));
        eH=max(abs(X(G.N+1:end)-yearRef(G.N+1:end)));
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
a=state(X(1:G.N),X(G.N+1:end),G,p);
out=struct('version','TH_v9_20260923_annular','T',X(1:G.N),'h',X(G.N+1:end), ...
    'Psoil',p.rho*p.g*X(G.N+1:end),'soilNodes',G.ids,'state',a,'pAll',G.xy,'mesh',S, ...
    'timeDay',t/p.day,'history',hist,'cycles',cycles,'snapshots',snaps,'parameters',p, ...
    'initialInfo',initInfo,'initialX',start,'periodicConverged',cycleConverged);
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
    'minimumGeometricEdgeWeight',G.minimumGeometricEdgeWeight);
out.leakage=struct('segmentRate_m2_per_s',leakHistSegment, ...
    'jointRate_m2_per_s',leakHistJoint,'bodyRate_m2_per_s',leakHistBody, ...
    'jointAnglesRad',p.leakage.jointAnglesRad(:)', ...
    'segmentCount',p.leakage.segmentCount,'interiorHead_m',p.leakage.interiorHead);
out.leakage=summarize_segment_leakage(out.history,out.leakage,t0/p.day);
if p.mechanics.enable, out.mechanics=mechanics(X,start,G,S,p); end
if strcmp(p.initial.mode,'spinup') && ~cycleConverged
    warning('TH:spinup','Maximum spin-up duration reached WITHOUT meeting periodic tolerances. Do not label the state independent of initialization.');
end
if p.saveOutput
    save([p.outputPrefix '_results.mat'],'out','-v7');
    checkpoint=makeCheckpoint(X,t,G,p,initInfo); %#ok<NASGU>
    save([p.outputPrefix '_checkpoint.mat'],'checkpoint','hist','cycles','snaps', ...
        'leakHistSegment','leakHistJoint','leakHistBody');
end
if p.plot, plotOutput(out,G); end
end

function p=defaults()
p.mesh=[]; p.meshFile=fullfile(fileparts(mfilename('fullpath')),'example_mesh_repaired_interface.mat');
p.day=86400; p.yearDays=365; p.durationDays=730;
p.dtInitial=.1*p.day; p.dtMax=p.day; p.dtMin=p.day/4096;
p.outputEveryDays=30; p.outputPrefix='TH_v9'; p.saveOutput=true; p.plot=true; p.verbose=true;
p.rho=1000; p.g=9.81; p.L=3.34e5; p.Tm=273.15;
p.porosity=.35; p.thetaR=.05; p.vgAlpha=.65; p.vgN=1.67; % alpha in 1/m, illustrative silt
p.Ksat=1e-14/1e-3*p.rho*p.g; % m/s, inherited intrinsic permeability
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
    'jointAnglesRad',[],'jointWidth_m',[],'jointConductance',[],'jointEntryHead',[]);
p.interfaceHTC=5;
p.advection=true;
p.initial=struct('mode','profile','surfaceT',-1,'gradient',.006,'T',[],'h',[], ...
    'checkpointFile','','startDay',0);
p.spinup=struct('tolT',.02,'tolH',.02,'tolWater',1e-4,'tolIce',1e-4, ...
    'consecutiveYears',2,'productionYears',1);
p.solver=struct('maxIter',40,'tol',1e-8,'stepT',2,'stepH',50,'fd',1e-8, ...
    'lineSearch',16,'massScale',.01,'energyScale',1e7, ...
    'jacobianReuse',4,'jacobianRefreshRatio',.65, ...
    'dtTargetResidualEvals',180,'dtGrowthSafety',1.10, ...
    'dtGrowthMax',1.8,'dtShrinkMin',.5);
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
assert(p.thetaR>=0 && p.thetaR<p.porosity && p.porosity<1);
assert(p.vgN>1 && p.vgAlpha>0 && p.Ksat>0 && p.Ss>0);
assert(p.dtMin>0 && p.dtInitial>=p.dtMin && p.dtMax>=p.dtInitial);
assert(p.outputEveryDays>0 && p.durationDays>0);
assert(p.topHTC>0 && p.tunnelHTC>0 && p.interfaceHTC>0);
assert(p.initial.gradient>=0 && p.rho>0 && p.L>0);
assert(any(strcmp(p.initial.mode,{'profile','fields','checkpoint','spinup'})));
assert(p.liningThickness>0 && p.liningK>=0 && p.iceImpedance>=0);
assert(p.solver.jacobianReuse>=0 && p.solver.jacobianReuse==fix(p.solver.jacobianReuse));
assert(p.solver.jacobianRefreshRatio>0 && p.solver.jacobianRefreshRatio<1);
assert(p.solver.dtTargetResidualEvals>0 && p.solver.dtGrowthSafety>0);
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
G.Kconc=sparse(G.N,G.N); G.vconc=zeros(G.N,1); G.qsets={}; G.ctri=zeros(0,3);
for name={'qTunnelGlobal','qInvertGlobal'}
    if isfield(S,name{1}) && ~isempty(S.(name{1}))
        qs=S.(name{1}); G.qsets{end+1}=qs;
        for e=1:size(qs,1)
            ids=qs(e,:); [K,v]=quadThermal(G.xy(ids,:));
            G.Kconc(ids,ids)=G.Kconc(ids,ids)+p.kc*K;
            G.vconc(ids)=G.vconc(ids)+v;
        end
    end
end
if isfield(S,'tInvertTriGlobal') && ~isempty(S.tInvertTriGlobal)
    G.ctri=S.tInvertTriGlobal; [ar,gr]=triGeometry(G.xy,G.ctri);
    [ii,jj,vv]=triTriplets(G.ctri,ar,gr);
    G.Kconc=G.Kconc+sparse(ii,jj,p.kc*vv,G.N,G.N);
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

function [K,v]=quadThermal(x)
K=zeros(4); v=zeros(4,1); signs=[];
for xi=[-1 1]/sqrt(3)
    for eta=[-1 1]/sqrt(3)
        [N,d]=qshape(xi,eta); J=d*x; dj=det(J); signs(end+1)=sign(dj); %#ok<AGROW>
        assert(abs(dj)>1e-12,'Degenerate Q4.'); dn=(J\d)';
        K=K+dn*dn'*abs(dj); v=v+N'*abs(dj);
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
% Split each interface edge at sector boundaries. Nodal binning alone gives
% a roughly 6% artificial sector bias on the supplied 51-edge circular mesh.
B.sectorNodeWeight=sparse(G.ns,p.leakage.segmentCount);
sectorAngle=2*pi/p.leakage.segmentCount;
for edge=1:size(G.interface,1)
    nodes=G.interface(edge,1:2); local=G.map(nodes);
    theta=mod(atan2(x(nodes,2)-G.center(2),x(nodes,1)-G.center(1)) ...
        -p.leakage.segmentStartAngleRad,2*pi);
    if theta(2)-theta(1)>pi, theta(2)=theta(2)-2*pi; end
    if theta(2)-theta(1)<-pi, theta(2)=theta(2)+2*pi; end
    assert(abs(diff(theta))>1e-12,'Degenerate angular interface edge.');
    low=min(theta); high=max(theta);
    crossings=(floor(low/sectorAngle)+1:ceil(high/sectorAngle)-1)*sectorAngle;
    fractions=sort([0;(crossings(:)-theta(1))/diff(theta);1]);
    len=norm(x(nodes(2),:)-x(nodes(1),:));
    for piece=1:numel(fractions)-1
        u=fractions(piece); v=fractions(piece+1); mid=(u+v)/2;
        sector=floor(mod(theta(1)+(u+v)/2*diff(theta),2*pi)/sectorAngle)+1;
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
            crossings=(floor((low-p.leakage.segmentStartAngleRad)/sectorAngle)+1: ...
                ceil((high-p.leakage.segmentStartAngleRad)/sectorAngle)-1)*sectorAngle ...
                +p.leakage.segmentStartAngleRad;
            cuts=[low,crossings,high];
            for piece=1:numel(cuts)-1
                fractions=sort(([cuts(piece),cuts(piece+1)]-theta(1))/diff(theta));
                midpoint=mean(fractions);
                len=norm(x(nodes(2),:)-x(nodes(1),:))*diff(fractions);
                sector=floor(mod(mean(cuts(piece:piece+1))- ...
                    p.leakage.segmentStartAngleRad,2*pi)/sectorAngle)+1;
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

function [X,t,info]=initialFields(G,B,p)
t=p.initial.startDay*p.day;
T=p.initial.surfaceT+p.initial.gradient*(G.xy(:,2)-G.z0);
h=G.xy(G.ids,2)-G.z0-waterTableDepthAt(t,p);
info=struct('mode',p.initial.mode,'note','Initial thermal and hydraulic transients are retained, not erased.');
if strcmp(p.initial.mode,'checkpoint')
    C=load(p.initial.checkpointFile,'checkpoint'); c=C.checkpoint;
    assert(strcmp(c.version,'TH_v9_20260923_annular'), ...
        'Incompatible checkpoint version.');
    assert(isequal(c.pAll,G.xy) && isequal(c.soilNodes,G.ids) && isequal(c.tSoil,G.tri) ...
        && isequal(c.qsets,G.qsets) && isequal(c.ctri,G.ctri) && isequal(c.interface,G.interface),'Checkpoint mesh mismatch.');
    assert(isequal(c.signature,signature(p)),'Checkpoint constitutive parameters differ. Start a fields run explicitly to change material.');
    T=c.X(1:G.N); h=c.X(G.N+1:end); t=c.time;
elseif strcmp(p.initial.mode,'fields')
    assert(numel(p.initial.T)==G.N && numel(p.initial.h)==G.ns,'Supply initial.T (all nodes) and initial.h (soil-node order).');
    T=p.initial.T(:); h=p.initial.h(:);
elseif ~isempty(p.initial.T) || ~isempty(p.initial.h)
    error('Use initial.mode=fields for supplied fields.');
end
assert(all(isfinite(T)) && all(T>-100 & T<100) && all(isfinite(h)));
X=[T;h];
% No resetting measured fields to BC values here: changed forcing starts at first step.
info.initialFarBoundaryMismatch=0;
if ~isempty(B.far)
    boundaryHead=farHeadAt(t,G,B,p);
    info.initialFarBoundaryMismatch=max(abs(h(B.far)-boundaryHead));
end
end

function s=signature(p)
% Include every physical/boundary setting that affects the continuation.
nonphysical={'mesh','meshFile','durationDays','dtInitial','dtMax','dtMin', ...
    'outputEveryDays','outputPrefix','saveOutput','plot','verbose','initial', ...
    'solver','spinup','mechanics'};
s=rmfield(p,nonphysical);
end

function c=makeCheckpoint(X,t,G,p,info)
c=struct('version','TH_v9_20260923_annular','X',X,'time',t,'pAll',G.xy,'soilNodes',G.ids, ...
    'tSoil',G.tri,'qsets',{G.qsets},'ctri',G.ctri,'interface',G.interface, ...
    'signature',signature(p),'parameters',p,'initialInfo',info);
end

function a=state(T,h,G,p)
% Painter-Karra beta=omega=1 limit, apparent liquid saturation closure (Eq 15).
% r=S*(-rho*g*h), F=S*(rho*L*max(-T,0)/Tm).
% sl=min(r,F); si=max(0,1-F/r). Exactly satisfies sl/(1-si)=r.
Ts=T(G.ids); r=retention(h,p);
hFreeze=-(p.L/(p.g*p.Tm))*max(-Ts,0); F=retention(hFreeze,p);
sl=min(r,F); si=max(0,1-F./r); sg=1-sl-si;
a.thetaW=p.porosity*sl; a.thetaI=p.porosity*si; a.thetaA=p.porosity*sg;
a.total=a.thetaW+a.thetaI;
a.elasticWater=p.Ss*max(h,0); % first-order compressibility storage, not pore occupancy
mass=a.total+a.elasticWater;
cs=(1-p.porosity)*p.rhoS*p.cs;
Es=cs*Ts+p.rho*(a.thetaW.*(p.L+p.cw*Ts)+a.thetaI.*p.ci.*Ts) ...
    +p.rho*a.elasticWater.*(p.L+p.cw*Ts);
a.E=G.vconc.*(p.rhoC*p.cc*T); a.E(G.ids)=a.E(G.ids)+G.vol.*Es;
a.W=G.vol.*mass;
a.k=(1-p.porosity)*p.ks+a.thetaW*p.kw+a.thetaI*p.ki+a.thetaA*p.ka;
% Mualem permeability based on liquid effective saturation plus ice impedance.
Se=(sl-p.thetaR/p.porosity)/(1-p.thetaR/p.porosity); Se=max(0,min(1,Se));
m=1-1/p.vgN; kr=sqrt(Se).*(1-(1-Se.^(1/m)).^m).^2;
a.K=p.Ksat*kr.*10.^(-p.iceImpedance*si);
a.sl=sl; a.si=si; a.sg=sg;
end

function s=retention(h,p)
m=1-1/p.vgN; sr=p.thetaR/p.porosity;
s=sr+(1-sr)*(1+(p.vgAlpha*max(-h,0)).^p.vgN).^(-m);
end

function [R,D]=residual(X,old,t,dt,G,B,p,scaled)
T=X(1:G.N); h=X(G.N+1:end); a=state(T,h,G,p);
a0=state(old(1:G.N),old(G.N+1:end),G,p);
ke=mean(reshape(a.k(G.local),size(G.local)),2); kt=sparse(G.ri,G.ci,G.kvals.*repmat(ke,9,1),G.N,G.N)+G.C;
Kh=mean(reshape(a.K(G.local),size(G.local)),2);
trans=accumarray(G.edgeMap,G.edgeBase.*Kh(G.edgeElement),[size(G.edges,1) 1]);
H=h-G.xy(G.ids,2); i=G.edges(:,1); j=G.edges(:,2); F=trans.*(H(i)-H(j));
qLeak=zeros(G.ns,1); qTop=qLeak;
if nargout>1
    leak=struct('body',0,'bySegment',zeros(p.leakage.segmentCount,1), ...
        'byJoint',zeros(numel(p.leakage.jointAnglesRad),1));
end
if p.tunnelDrain
    if nargout>1, [qLeak,leak]=liningLeak(h,a,B,p);
    else, qLeak=liningLeak(h,a,B,p); end
end
if p.topSeepage
    qTop=p.topLeakConductance*B.topW(G.ids).*max(h,0);
end
rm=(a.W-a0.W)/dt+G.inc*F+qLeak+qTop;
qFar=zeros(G.ns,1); qFar(B.far)=-rm(B.far); % exact prescribed-head reaction, outward positive
re=(a.E-a0.E)/dt+kt*T;
if p.advection
    tw=T(G.ids); upstream=tw(j); upstream(F>=0)=tw(i(F>=0));
    Fa=p.rho*(p.L+p.cw*upstream).*F;
    re(G.ids)=re(G.ids)+G.inc*Fa;
    tfar=p.groundwaterInflowT+p.groundwaterInflowGradient*(G.xy(G.ids,2)-G.z0);
    Tadv=tw; Tadv(qFar<0)=tfar(qFar<0);
    re(G.ids)=re(G.ids)+p.rho*(p.L+p.cw*tw).*(qLeak+qTop) ...
        +p.rho*(p.L+p.cw*Tadv).*qFar;
end
[Ttop,Ttun]=thermalValues(t,p);
qHeat=-p.geothermalFlux*B.bottomW;
if isfinite(p.topHTC), qHeat=qHeat+p.topHTC*B.topW.*(T-Ttop); end
if isfinite(p.tunnelHTC), qHeat=qHeat+p.tunnelHTC*B.tunnelW.*(T-Ttun); end
re=re+qHeat;
qTdir=zeros(G.N,1); qTdir(B.Tdir)=-re(B.Tdir);
if nargout>1
    D.water=sum(a.W); D.energy=sum(a.E);
    D.Qout=sum(qFar+qLeak+qTop); D.leak=sum(qLeak);
    D.leakBody=leak.body; D.leakBySegment=leak.bySegment;
    D.leakByJoint=leak.byJoint;
    D.Eout=sum(qHeat+qTdir);
    if p.advection
        D.Eout=D.Eout+sum(p.rho*(p.L+p.cw*tw).*(qLeak+qTop)+p.rho*(p.L+p.cw*Tadv).*qFar);
    end
    D.waterError=sum(a.W-a0.W)+dt*D.Qout;
    D.energyError=sum(a.E-a0.E)+dt*D.Eout;
    D.minTrans=min(trans); D.waterReaction=qFar; D.F=F;
end
R=dt*[re;rm];
if scaled, R=R./G.rs; end
% Boundary rows always dimensionless, same in numerical Jacobian and line search.
R(G.N+B.far)=(h(B.far)-farHeadAt(t,G,B,p))/10;
if isinf(p.topHTC), R(B.top)=T(B.top)-Ttop; end
if isinf(p.tunnelHTC), ids=find(B.tunnelW>0); R(ids)=T(ids)-Ttun; end
end

function [qLeak,ledger]=liningLeak(h,~,B,p)
% Annular body Robin conductance and calibrated finite-width joint bands.
% Each band flow is distributed to adjacent soil nodes, preserving the
% total water and energy budgets on this 2D slice.
% The driving potential crosses the real lining thickness in gravity:
% H_out-H_in=(h_out-z_out)-(h_in-z_in).
drivingHead=h-p.leakage.interiorHead+B.liningZDrop;
bodyFlux=p.liningK*B.bodyRadialFactor.*max(drivingHead,0);
body=B.bodyW.*bodyFlux;
qLeak=body;
nj=numel(p.leakage.jointAnglesRad);
joint=zeros(nj,1);
if nj>0
    c=p.leakage.jointConductance(:);
    if isscalar(c), c=repmat(c,nj,1); end
    entry=p.leakage.jointEntryHead(:);
    if isempty(entry), entry=zeros(nj,1);
    elseif isscalar(entry), entry=repmat(entry,nj,1); end
    nodeFlow=c(B.jointID).*B.jointNormalizedWeight.* ...
        max(drivingHead(B.jointNodeIndex)-entry(B.jointID),0);
    joint=accumarray(B.jointID,nodeFlow,[nj 1]);
    qLeak=qLeak+accumarray(B.jointNodeIndex,nodeFlow,[numel(h) 1]);
end
if nargout>1
    bySegment=B.bodySectorNodeWeight'*bodyFlux;
    if nj>0
        for sector=1:p.leakage.segmentCount
            jointID=B.jointSectorID{sector};
            nodeID=B.jointSectorNodeIndex{sector};
            weight=B.jointSectorNormalizedWeight{sector};
            bySegment(sector)=bySegment(sector)+sum(c(jointID).*weight.* ...
                max(drivingHead(nodeID)-entry(jointID),0));
        end
    end
    ledger=struct('body',sum(body),'bySegment',bySegment,'byJoint',joint);
end
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

function [X,info]=step(old,t,dt,G,B,p)
% Adaptive modified Newton (chord Newton).  A sparse coloured-FD tangent is
% reused while it remains predictive.  Poor contraction forces a rebuild;
% a failed line search with a reused tangent gets one immediate fresh-tangent
% retry.  Therefore reuse changes cost, not the converged discrete equations.
X=old; J=[]; F=[]; Jage=inf;
info=struct('ok',false,'iter',0,'norm',inf,'reason','iteration limit', ...
    'jacobianBuilds',0,'linearSolves',0,'residualEvaluations',0);
for k=1:p.solver.maxIter
    R=residual(X,old,t+dt,dt,G,B,p,true); nr=norm(R,inf);
    info.residualEvaluations=info.residualEvaluations+1;
    info.iter=k; info.norm=nr;
    if ~isfinite(nr), info.reason='nonfinite residual'; return; end
    if nr<p.solver.tol, info.ok=true; info.reason='converged'; return; end
    accepted=false;
    if isempty(J) || Jage>=p.solver.jacobianReuse
        [J,nfd]=jacobian(X,old,t+dt,dt,G,B,p,R);
        F=factorize(J);
        info.jacobianBuilds=info.jacobianBuilds+1;
        info.residualEvaluations=info.residualEvaluations+nfd; Jage=0;
    end
    merit0=norm(R,2);
    for tangentAttempt=1:2
        dy=-(F\R); info.linearSolves=info.linearSolves+1; dx=G.xs.*dy;
        if any(~isfinite(dx)), info.reason='singular tangent'; return; end
        lam=min([1,p.solver.stepT/max(max(abs(dx(1:G.N))),eps), ...
            p.solver.stepH/max(max(abs(dx(G.N+1:end))),eps)]);
        for ls=1:p.solver.lineSearch
            candidate=X+lam*dx;
            if all(candidate(1:G.N)>-100 & candidate(1:G.N)<100)
                rr=residual(candidate,old,t+dt,dt,G,B,p,true);
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
            [J,nfd]=jacobian(X,old,t+dt,dt,G,B,p,R);
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

function [J,nEval]=jacobian(X,old,t,dt,G,B,p,R)
v=zeros(numel(G.jrow),1);
for c=1:numel(G.groups)
    cols=G.groups{c}; delta=p.solver.fd*max(abs(X(cols)),G.xs(cols));
    xp=X; xp(cols)=xp(cols)+delta;
    dr=residual(xp,old,t,dt,G,B,p,true)-R;
    nz=G.jnz{c}; ii=G.jrow(nz); jj=G.jcol(nz);
    v(nz)=dr(ii)./delta(G.jloc{c}).*G.xs(jj);
end
J=sparse(G.jrow,G.jcol,v,numel(X),numel(X));
nEval=numel(G.groups);
end

function q=darcy(X,G,p)
a=state(X(1:G.N),X(G.N+1:end),G,p); h=X(G.N+1:end);
gradH=zeros(numel(G.area),2);
for j=1:3, gradH=gradH+G.grad(:,:,j).*h(G.local(:,j)); end
gradH(:,2)=gradH(:,2)-1;
q.elemFlux=-mean(reshape(a.K(G.local),size(G.local)),2).*gradH;
q.elemCentroid=(G.xy(G.tri(:,1),:)+G.xy(G.tri(:,2),:)+G.xy(G.tri(:,3),:))/3;
q.head=h-G.xy(G.ids,2); q.note='Darcy centroid flux; budgets use exact assembled edge/boundary fluxes.';
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
p=defaults(); p.plot=false; p.saveOutput=false; p.verbose=false;
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
p=defaults(); p.mechanics.enable=false;
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
report=struct('passed',true,'jacobianRelativeError',jacerr,'sealedWaterDrift_m2',sum(a.W)-W0, ...
    'maxStepWaterDefect_m2',maxW,'maxStepEnergyDefect_J_per_m',maxE,'sealedMaxIce',max(a.thetaI), ...
    'sectorLengthSpread_m',max(segmentLengths)-min(segmentLengths), ...
    'finiteWidthJointFlow_m2_per_s',ledger.byJoint);
disp(report);
end
