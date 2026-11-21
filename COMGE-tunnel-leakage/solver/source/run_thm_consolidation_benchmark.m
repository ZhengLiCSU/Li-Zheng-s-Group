function study=run_thm_consolidation_benchmark(options)
%RUN_THM_CONSOLIDATION_BENCHMARK Saturated isothermal Biot sine-mode patch.
% Rectangular column [0,W] x [0,H], y=z positive down. All ux are fixed;
% bottom uz is fixed; the top has zero incremental TOTAL normal traction.
% Hydraulic top excess head is zero, sides/bottom have zero external flow.
% hAbsolute=z+hBaseline+hExcess stays positive. The baseline pressure and
% selfweight are included in an installed initial equilibrium, then cancel
% from the incremental mechanical equations. No phase switches are active.
%
% Initial excess head=A*sin(pi*z/(2H)). With M=E(1-nu)/((1+nu)(1-2nu)),
%   S_eff=Ss+rho*g*alphaB^2/M, c=Ksat/S_eff,
%   hExcess(z,t)=A*exp[-c*(pi/(2H))^2*t]*sin(pi*z/(2H)),
%   uzIncrement(z,t)=alphaB*rho*g*A/M * (1-decay)*cos(pi*z/(2H))/(pi/(2H)).
% Thus positive pressure loss gives downward settlement at the free top.
%
% The solved unknown is [hExcess;uIncrement]. Water storage uses the same
% median-SCV owner deposit, P1 pressure interpolation and element div(u) as
% v16: M=phi0+alphaB*div(uIncrement)+Ss*hAbsolute in this saturated limit.
% Mechanics uses soil_poroelastic_fem; flux uses soil_cvfe_flux. One coupled
% sparse block is solved at every backward-Euler step, with no staggered
% pressure/displacement update. The SCV mass and pressure-load blocks are
% verified as adjoints (Bpressure=rho*g*Mdisplacement').
%
% This verifies the PROJECT'S COUPLING KERNELS in the linear warm saturated
% limit, not the complete tunnel v16 solver, freezing, contact or joint law.
% Default H=20 m retains a measurable mode after 1 day (rather than complete
% diffusion at H=1 m). Vertical node counts=[5 9 17], 3 horizontal nodes,
% dt=10 s; the finest case repeats at 5 s to measure temporal contamination.
% Spatial L2 is normalized by rectangle area (RMS head error [m]) and uses
% independent degree-5 triangle quadrature of the continuous exact field.
%
% Options are the fields of configuration below. No MATLAB toolboxes needed.
if nargin<1||isempty(options), options=struct(); end
root=fileparts(mfilename('fullpath'));
c=struct('verticalNodes',[5 9 17],'horizontalNodes',3, ...
    'height_m',20,'width_m',2,'headAmplitude_m',.1,'baselineHead_m',2, ...
    'durationDays',1,'dtSeconds',10,'checkTimeRefinement',true, ...
    'E_Pa',80e6,'poissonRatio',.35,'biot',.65, ...
    'Ksat_m_s',1e-7,'specificStorage_per_m',1e-5, ...
    'saveArtifacts',true,'verbose',true, ...
    'outputFolder',fullfile(root,'thm_consolidation_benchmark'));
assert(isstruct(options)&&isscalar(options),'THMConsolidation:Options', ...
    'Supply a scalar options structure.');
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(c,names{k}),'THMConsolidation:Options', ...
        'Unknown benchmark option: %s.',names{k});
    c.(names{k})=options.(names{k});
end
validateConfiguration(c); c.verticalNodes=c.verticalNodes(:)';
rho=1000; gravity=9.81;
M=c.E_Pa*(1-c.poissonRatio)/((1+c.poissonRatio)*(1-2*c.poissonRatio));
coupledStorage=rho*gravity*c.biot^2/M;
effectiveStorage=c.specificStorage_per_m+coupledStorage;
theory=struct('constrainedModulus_Pa',M, ...
    'specificStorage_per_m',c.specificStorage_per_m, ...
    'poroelasticStorage_per_m',coupledStorage, ...
    'effectiveStorage_per_m',effectiveStorage, ...
    'diffusivity_m2_s',c.Ksat_m_s/effectiveStorage, ...
    'waveNumber_per_m',pi/(2*c.height_m),'rho_kg_m3',rho,'gravity_m_s2',gravity);
results=cell(numel(c.verticalNodes),1); metrics=struct([]); started=tic;
for k=1:numel(c.verticalNodes)
    results{k}=runCase(c.verticalNodes(k),c.dtSeconds,c,theory);
    if k==1, metrics=results{k}.metrics; else, metrics(k,1)=results{k}.metrics; end %#ok<AGROW>
    if c.verbose
        r=results{k}.metrics;
        fprintf('Consolidation nz=%d, dt=%g s: RMS head %.6g m, budget %.3g, %.3g s.\n', ...
            r.verticalNodes,r.dtSeconds,r.rmsHeadError_m,r.relativeCumulativeWaterDefect,r.wallSeconds);
    end
end
comparison=struct2table(metrics);
comparison.observedOrder=[NaN;log(comparison.rmsHeadError_m(1:end-1)./ ...
    comparison.rmsHeadError_m(2:end))./log(comparison.dz_m(1:end-1)./comparison.dz_m(2:end))];
timeRefinement=struct('performed',false);
if c.checkTimeRefinement
    fine=runCase(c.verticalNodes(end),c.dtSeconds/2,c,theory);
    coarse=results{end};
    difference=fieldDifference(coarse.headExcess_m,fine.headExcess_m,fine.xy,fine.tri);
    timeRefinement=struct('performed',true,'fine',fine, ...
        'rmsHeadChange_m',difference,'estimatedFullStepTemporalError_m',2*difference, ...
        'estimatedTemporalFractionOfFineL2',2*difference/max(fine.metrics.rmsHeadError_m,realmin), ...
        'note','Backward-Euler dt versus dt/2 estimate at fixed mesh; numerical evidence, not a rigorous time-error bound.');
    if c.verbose
        fprintf('Consolidation dt/2: RMS change %.3g m, estimated temporal/spatial fraction %.3g.\n', ...
            difference,timeRefinement.estimatedTemporalFractionOfFineL2);
    end
end
study=struct('benchmark','Saturated isothermal plane-strain Biot consolidation sine mode', ...
    'version','THM_consolidation_kernel_20260930','configuration',c, ...
    'theory',theory,'caseResults',{results},'comparison',comparison, ...
    'timeRefinement',timeRefinement,'totalWallSeconds',toc(started), ...
    'scope',['Monolithic linear warm-saturated kernel verification using project SCV storage, ' ...
        'physical Darcy subfaces and soil elasticity. It does not run the complete tunnel v16 ' ...
        'solver or validate frozen soil, contact, density expansion, gasket or field parameters.']);
if c.saveArtifacts
    folder=char(c.outputFolder); if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,'thm_consolidation_results.mat'),'study','-v7');
    writetable(comparison,fullfile(folder,'thm_consolidation_mesh_errors.csv'));
    if timeRefinement.performed
        temporal=struct2table(timeRefinement.fine.metrics);
        temporal.rmsChangeFromFullStep_m=timeRefinement.rmsHeadChange_m;
        temporal.estimatedTemporalFraction=timeRefinement.estimatedTemporalFractionOfFineL2;
        writetable(temporal,fullfile(folder,'thm_consolidation_time_check.csv'));
    end
    writeReport(study,fullfile(folder,'thm_consolidation_report.txt'));
end
end

function result=runCase(nz,requestedDt,c,theory)
started=tic; nx=c.horizontalNodes;
[x,z]=ndgrid(linspace(0,c.width_m,nx),linspace(0,c.height_m,nz));
xy=[x(:),z(:)]; ns=size(xy,1); tri=zeros(2*(nx-1)*(nz-1),3); cursor=0;
node=reshape(1:ns,nx,nz);
for j=1:nz-1
    for i=1:nx-1
        a=node(i,j); b=node(i+1,j); d=node(i,j+1); e=node(i+1,j+1);
        if mod(i+j,2)==0, pair=[a b d;b e d]; else, pair=[a b e;a e d]; end
        tri(cursor+(1:2),:)=pair; cursor=cursor+2;
    end
end
ne=size(tri,1); top=find(xy(:,2)==0); bottom=find(xy(:,2)==c.height_m);
initialHead=c.headAmplitude_m*sin(theory.waveNumber_per_m*xy(:,2));
initialHead(top)=0;
baselineHead=xy(:,2)+c.baselineHead_m;
p=struct('rho',theory.rho_kg_m3,'g',theory.gravity_m_s2,'L',3.35e5,'Tm',273.15, ...
    'porosity',.35,'thetaR',.05,'vgAlpha',.65,'vgN',1.67, ...
    'Ksat',c.Ksat_m_s,'Ss',c.specificStorage_per_m,'iceImpedance',0, ...
    'cw',4180,'ci',2100,'rhoS',2650,'cs',800,'kw',.57,'ki',2.2,'ks',2.3,'ka',.025);
settings=struct('E',c.E_Pa,'nu',c.poissonRatio,'biot',c.biot,'alphaThermal',0, ...
    'rhoBulk',(1-p.porosity)*p.rhoS+p.porosity*p.rho, ...
    'gravity',[0;p.g],'phasePressure','liquid','phaseParameters',p);
mechanics=soil_poroelastic_fem(xy,tri,5,baselineHead+initialHead,5,settings);
Q=build_soil_scv_quadrature(xy,tri,5,0);
cvfe=build_soil_cvfe_flux(xy,tri,(1:ns)');
[~,fluxH]=soil_cvfe_flux(baselineHead,p.Ksat*ones(ne,1),cvfe,'element');
hydraulic=cvfe.inc*fluxH;
elementSelect=sparse((1:numel(Q.element))',Q.element,1,numel(Q.element),ne);
pressureStorage=p.Ss*(Q.deposit*Q.P);
displacementStorage=c.biot*(Q.deposit*elementSelect*mechanics.V);
pressureLoad=p.rho*p.g*mechanics.A_pressure*(Q.elementAverage*Q.P);
adjointError=norm(pressureLoad-p.rho*p.g*displacementStorage',inf) ...
    /max(norm(pressureLoad,inf),realmin);
fixed=unique([2*(1:ns)'-1;2*bottom]);
free=setdiff((1:2*ns)',fixed);
referenceU=zeros(2*ns,1);
referenceForce=mechanics.fbody+mechanics.fpore;
referenceU(free)=mechanics.K(free,free)\referenceForce(free);
referenceResidual=mechanics.K*referenceU-referenceForce;
referenceRelativeResidual=norm(referenceResidual(free),inf)/max(1,norm(referenceForce(free),inf));
endSeconds=c.durationDays*86400; steps=ceil(endSeconds/requestedDt); dt=endSeconds/steps;
system=[pressureStorage+dt*hydraulic,displacementStorage; ...
    -pressureLoad,mechanics.K];
fixedRows=[top;ns+fixed];
system(fixedRows,:)=0;
system=system+sparse(fixedRows,fixedRows,1,3*ns,3*ns);
% Row and column scaling change conditioning, not the coupled equations.
uScale=max(c.biot*p.rho*p.g/theory.constrainedModulus_Pa ...
    *c.headAmplitude_m*c.height_m,1e-9);
columnScale=[c.headAmplitude_m*ones(ns,1);uScale*ones(2*ns,1)];
scaled=system*spdiags(columnScale,0,3*ns,3*ns);
rowScale=1./full(max(abs(scaled),[],2));
scaled=spdiags(rowScale,0,3*ns,3*ns)*scaled;
factor=decomposition(scaled,'lu');
head=initialHead; u=zeros(2*ns,1);
initialMass=Q.deposit*(p.porosity+p.Ss*(Q.P*(baselineHead+initialHead)));
cumulativeDefect=0; absoluteDefect=0; exportedWater=0; absoluteExport=0;
maxLocalDefect=0; maxMechanicalResidual=0; minimumHead=min(baselineHead+head);
maxVolumeStrain=0; history=zeros(steps,5);
mechanicalRight=-pressureLoad*initialHead;
forceScale=max(1,norm(pressureLoad*initialHead,inf));
for step=1:steps
    oldHead=head; oldU=u;
    rhs=[pressureStorage*head+displacementStorage*u;mechanicalRight];
    rhs(fixedRows)=0;
    solution=columnScale.*(factor\(rowScale.*rhs));
    head=solution(1:ns); u=solution(ns+1:end);
    % Evaluate the same physical flux kernel, not an inferred boundary flux.
    flow=soil_cvfe_flux(baselineHead+head,p.Ksat*ones(ne,1),cvfe,'element');
    waterIncrement=pressureStorage*(head-oldHead)+displacementStorage*(u-oldU);
    waterResidual=waterIncrement+dt*(cvfe.inc*flow);
    drainReaction=zeros(ns,1); drainReaction(top)=-waterResidual(top)/dt;
    localDefect=waterResidual+dt*drainReaction;
    budget=sum(waterIncrement)+dt*sum(drainReaction);
    cumulativeDefect=cumulativeDefect+budget; absoluteDefect=absoluteDefect+abs(budget);
    exportedWater=exportedWater+dt*sum(drainReaction);
    absoluteExport=absoluteExport+dt*sum(abs(drainReaction));
    maxLocalDefect=max(maxLocalDefect,max(abs(localDefect)));
    mechanicalResidual=mechanics.K*u-pressureLoad*(head-initialHead);
    maxMechanicalResidual=max(maxMechanicalResidual,norm(mechanicalResidual(free),inf));
    minimumHead=min(minimumHead,min(baselineHead+head));
    maxVolumeStrain=max(maxVolumeStrain,max(abs(mechanics.V*u)));
    history(step,:)=[step*dt,sum(drainReaction),sum(waterIncrement),budget, ...
        norm(mechanicalResidual(free),inf)];
end
decay=exp(-theory.diffusivity_m2_s*theory.waveNumber_per_m^2*endSeconds);
exactHead=c.headAmplitude_m*decay*sin(theory.waveNumber_per_m*xy(:,2));
exactVertical=c.biot*p.rho*p.g/theory.constrainedModulus_Pa*c.headAmplitude_m ...
    *(1-decay)/theory.waveNumber_per_m*cos(theory.waveNumber_per_m*xy(:,2));
errorHead=exactError(head,xy,tri,@(z)c.headAmplitude_m*decay*sin(theory.waveNumber_per_m*z));
errorU=exactError(u(2:2:end),xy,tri,@(z)c.biot*p.rho*p.g/theory.constrainedModulus_Pa ...
    *c.headAmplitude_m*(1-decay)/theory.waveNumber_per_m*cos(theory.waveNumber_per_m*z));
initialError=exactError(initialHead,xy,tri,@(z)c.headAmplitude_m*sin(theory.waveNumber_per_m*z));
initialMassChange=sum(pressureStorage*(head-initialHead)+displacementStorage*u);
finalMass=initialMass+pressureStorage*(head-initialHead)+displacementStorage*u;
metrics=struct('verticalNodes',nz,'horizontalNodes',nx,'soilNodes',ns,'elements',ne, ...
    'dz_m',c.height_m/(nz-1),'dtSeconds',dt,'finalDays',endSeconds/86400, ...
    'acceptedSteps',steps,'rmsHeadError_m',errorHead,'nodalHeadLinf_m',max(abs(head-exactHead)), ...
    'relativeL2HeadError',errorHead/max(c.headAmplitude_m*decay/sqrt(2),realmin), ...
    'initialInterpolationRMS_m',initialError,'rmsVerticalIncrementError_m',errorU, ...
    'exactDecayFactor',decay,'exactTopSettlement_m',exactVertical(top(1)), ...
    'meanComputedTopSettlement_m',mean(u(2*top)), ...
    'minimumAbsolutePressureHead_m',minimumHead,'maxAbsIncrementalVolumeStrain',maxVolumeStrain, ...
    'maxTopExcessHead_m',max(abs(head(top))), ...
    'storageLoadAdjointRelativeError',adjointError, ...
    'installedReferenceRelativeResidual',referenceRelativeResidual, ...
    'maxFreeMechanicalResidual_N_per_m',maxMechanicalResidual, ...
    'relativeMechanicalResidual',maxMechanicalResidual/forceScale, ...
    'cumulativeWaterDefect_m2',cumulativeDefect,'sumAbsWaterDefect_m2',absoluteDefect, ...
    'maxLocalWaterDefect_m2',maxLocalDefect,'totalDrainedWater_m2',exportedWater, ...
    'totalStorageChange_m2',initialMassChange, ...
    'relativeCumulativeWaterDefect',abs(cumulativeDefect)/max(absoluteExport,realmin), ...
    'relativeSumAbsWaterDefect',absoluteDefect/max(absoluteExport,realmin), ...
    'wallSeconds',toc(started));
result=struct('xy',xy,'tri',tri,'headExcess_m',head,'absolutePressureHead_m',baselineHead+head, ...
    'initialHeadExcess_m',initialHead,'uIncrement_m',u,'installedReferenceU_m',referenceU, ...
    'initialWaterMass_m2',initialMass,'finalWaterMass_m2',finalMass, ...
    'exactNodalHead_m',exactHead,'exactNodalVerticalIncrement_m',exactVertical, ...
    'metrics',metrics,'history',history, ...
    'historyColumns',{{'time_s','drainOutflow_m2_s','waterIncrement_m2','waterDefect_m2','freeMechanicalResidual_N_m'}});
end

function value=exactError(field,xy,tri,exact)
[N,w]=triangleRule(); accumulated=0; totalArea=0;
for e=1:size(tri,1)
    ids=tri(e,:); x=xy(ids,:);
    area=abs(det([x(2,:)-x(1,:);x(3,:)-x(1,:)]))/2;
    q=N*x; difference=N*field(ids)-exact(q(:,2));
    accumulated=accumulated+area*sum(w.*difference.^2); totalArea=totalArea+area;
end
value=sqrt(accumulated/totalArea);
end

function value=fieldDifference(a,b,xy,tri)
[N,w]=triangleRule(); accumulated=0; totalArea=0;
for e=1:size(tri,1)
    ids=tri(e,:); x=xy(ids,:);
    area=abs(det([x(2,:)-x(1,:);x(3,:)-x(1,:)]))/2;
    difference=N*(a(ids)-b(ids));
    accumulated=accumulated+area*sum(w.*difference.^2); totalArea=totalArea+area;
end
value=sqrt(accumulated/totalArea);
end

function [N,w]=triangleRule()
a=.470142064105115; b=.059715871789770; c=.101286507323456; d=.797426985353087;
N=[1/3 1/3 1/3;a a b;a b a;b a a;c c d;c d c;d c c];
w=[.225;repmat(.132394152788506,3,1);repmat(.125939180544827,3,1)]; w=w/sum(w);
end

function validateConfiguration(c)
positive={'height_m','width_m','headAmplitude_m','baselineHead_m','durationDays', ...
    'dtSeconds','E_Pa','Ksat_m_s','specificStorage_per_m'};
for k=1:numel(positive)
    value=c.(positive{k});
    assert(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>0, ...
        'THMConsolidation:Options','%s must be a positive scalar.',positive{k});
end
assert(isnumeric(c.verticalNodes)&&isreal(c.verticalNodes)&&isvector(c.verticalNodes) ...
    &&numel(c.verticalNodes)>=2&&all(isfinite(c.verticalNodes)) ...
    &&all(c.verticalNodes>=3&c.verticalNodes==fix(c.verticalNodes)) ...
    &&all(diff(c.verticalNodes)>0),'THMConsolidation:Options', ...
    'verticalNodes must contain at least two increasing integer counts >=3.');
assert(isnumeric(c.horizontalNodes)&&isscalar(c.horizontalNodes)&&isfinite(c.horizontalNodes) ...
    &&c.horizontalNodes>=2&&c.horizontalNodes==fix(c.horizontalNodes), ...
    'THMConsolidation:Options','horizontalNodes must be an integer >=2.');
assert(isnumeric(c.poissonRatio)&&isreal(c.poissonRatio)&&isscalar(c.poissonRatio) ...
    &&isfinite(c.poissonRatio)&&c.poissonRatio>-1&&c.poissonRatio<.5 ...
    &&isnumeric(c.biot)&&isreal(c.biot)&&isscalar(c.biot)&&isfinite(c.biot) ...
    &&c.biot>=0&&c.biot<=1,'THMConsolidation:Options', ...
    'Require -1<nu<0.5 and 0<=biot<=1.');
for name={'checkTimeRefinement','saveArtifacts','verbose'}
    value=c.(name{1}); assert(isscalar(value)&&isreal(value)&&isfinite(value)&&ismember(value,[0 1]), ...
        'THMConsolidation:Options','Benchmark switches must be scalar logical 0/1.');
end
assert((ischar(c.outputFolder)&&isrow(c.outputFolder)&&~isempty(c.outputFolder)) ...
    ||(isstring(c.outputFolder)&&isscalar(c.outputFolder)&&strlength(c.outputFolder)>0), ...
    'THMConsolidation:Options','Supply a nonempty output-folder path.');
end

function writeReport(study,path)
fid=fopen(path,'w'); assert(fid>=0,'THMConsolidation:Report','Could not open report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
c=study.configuration; t=study.theory;
fprintf(fid,'%s\n\n',study.benchmark);
fprintf(fid,'Domain: width %.9g m, height %.9g m; z positive down.\n',c.width_m,c.height_m);
fprintf(fid,'ux fixed everywhere; uz bottom fixed; incremental total traction free at top.\n');
fprintf(fid,'Top excess head=0; bottom and sides natural zero external flow.\n');
fprintf(fid,'Initial h_excess=A sin(pi*z/(2H)); h_absolute=z+hBaseline+h_excess>0.\n');
fprintf(fid,'Initial gravity+pressure installation equilibrium is retained as the displacement reference.\n');
fprintf(fid,'Mcon=E(1-nu)/[(1+nu)(1-2nu)]=%.12g Pa.\n',t.constrainedModulus_Pa);
fprintf(fid,'S_eff=Ss+rho*g*alphaB^2/Mcon=%.12g 1/m; c=K/S_eff=%.12g m2/s.\n', ...
    t.effectiveStorage_per_m,t.diffusivity_m2_s);
fprintf(fid,'Exact head=A exp[-c(pi/(2H))^2 t] sin(pi*z/(2H)).\n');
fprintf(fid,'Exact downward settlement increment=alphaB*rho*g*A/Mcon*(1-decay)*cos(pi*z/(2H))/(pi/(2H)).\n');
fprintf(fid,'Monolithic [h_excess;u_increment], backward Euler, project T3 mechanics and median-SCV water storage/flux.\n');
fprintf(fid,'Water increments omit only algebraically cancelling constant baseline inventory.\n');
fprintf(fid,'Area-normalized L2/RMS evaluates the FE interpolant against the continuous exact field.\n\n');
fprintf(fid,'nz   dz(m)   dt(s)   RMShead(m)   order   relWaterBudget   wall(s)\n');
for k=1:height(study.comparison)
    r=study.comparison(k,:);
    fprintf(fid,'%d   %.8g   %.8g   %.10g   %.6g   %.4g   %.4g\n', ...
        r.verticalNodes,r.dz_m,r.dtSeconds,r.rmsHeadError_m,r.observedOrder, ...
        r.relativeCumulativeWaterDefect,r.wallSeconds);
end
if study.timeRefinement.performed
    r=study.timeRefinement;
    fprintf(fid,'\nFinest dt/2 RMS change=%.12g m; estimated full-step time/fine-space error fraction=%.6g.\n', ...
        r.rmsHeadChange_m,r.estimatedTemporalFractionOfFineL2);
end
fprintf(fid,'\nScope: %s\n',study.scope);
fprintf(fid,'CPU time is reported separately and is not an accuracy/pass criterion.\n');
end
