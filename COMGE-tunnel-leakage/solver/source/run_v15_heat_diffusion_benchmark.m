function study=run_v15_heat_diffusion_benchmark(options)
%RUN_V15_HEAT_DIFFUSION_BENCHMARK Analytic warm transient diffusion patch.
% Unit square, Z positive down: T(0,t)=2, other sides insulated. Initially
% T=2+sin(pi*z/2); exact T=2+exp(-alpha*(pi/2)^2*t)*sin(pi*z/2).
% All soil is warm and saturated. h=z+5 is hydrostatic, no tunnel or top
% seepage is present, and advection stays ENABLED to check spurious flow.
% This is verification of the discrete diffusion operator and storage,
% not validation of freeze/thaw constitutive parameters against experiment.
% Default: nx=[5 9 17] points per side, one day, fixed 120 s steps; the
% finest mesh is repeated at 60 s to quantify temporal contamination.
% Errors integrate the FE interpolant against the continuous exact field
% with independent positive degree-5 triangular quadrature.
% Optional fields: nx, durationDays, dtSeconds, checkTimeRefinement,
% saveArtifacts, outputFolder, verbose. Solver output is held in STUDY.
if nargin<1||isempty(options), options=struct(); end
root=fileparts(mfilename('fullpath'));
config=struct('nx',[5 9 17],'durationDays',1,'dtSeconds',120, ...
    'checkTimeRefinement',true,'saveArtifacts',true, ...
    'outputFolder',fullfile(root,'v15_heat_diffusion_benchmark'), 'verbose',true);
assert(isstruct(options)&&isscalar(options),'Supply a scalar options structure.');
names=fieldnames(options);
for j=1:numel(names)
    assert(isfield(config,names{j}),'Unknown benchmark option: %s.',names{j});
    config.(names{j})=options.(names{j});
end
assert(isnumeric(config.nx)&&isvector(config.nx)&&numel(config.nx)>=2&& ...
    all(isfinite(config.nx))&&all(config.nx>=3)&& ...
    all(config.nx==fix(config.nx))&&all(diff(config.nx)>0), ...
    'nx must contain at least two increasing integer point counts >=3.');
config.nx=config.nx(:)';
assert(isscalar(config.durationDays)&&isfinite(config.durationDays)&&config.durationDays>0);
assert(isscalar(config.dtSeconds)&&isfinite(config.dtSeconds)&&config.dtSeconds>0);
for flag={'checkTimeRefinement','saveArtifacts','verbose'}
    v=config.(flag{1}); assert(isscalar(v)&&any(v==[0 1]));
end
assert((ischar(config.outputFolder)&&isrow(config.outputFolder))|| ...
    (isstring(config.outputFolder)&&isscalar(config.outputFolder)));
caseResults=cell(numel(config.nx),1); rows=struct([]);
started=tic;
for j=1:numel(config.nx)
    caseResults{j}=runCase(config.nx(j),config.dtSeconds,config.durationDays);
    if j==1
        rows=caseResults{j}.metrics;
    else
        rows(j,1)=caseResults{j}.metrics; %#ok<AGROW>
    end
    if config.verbose
        r=caseResults{j}.metrics;
        fprintf('Warm diffusion nx=%d, dt=%g s: L2 %.6g C; max Darcy %.3g m/s; %.2g s.\n', ...
            r.nx,r.dtSeconds,r.L2_C,r.maxDarcyVelocity_m_s,r.wallSeconds);
    end
end
comparison=struct2table(rows);
comparison.observedOrder=[NaN;log(comparison.L2_C(1:end-1)./comparison.L2_C(2:end)) ./ ...
    log(comparison.dx_m(1:end-1)./comparison.dx_m(2:end))];
timeRefinement=struct('performed',false);
if config.checkTimeRefinement
    fine=runCase(config.nx(end),config.dtSeconds/2,config.durationDays);
    coarse=caseResults{end};
    change=finiteElementDifference(coarse.out.T,fine.out.T,fine.out.mesh);
    timeRefinement=struct('performed',true,'fine',fine, ...
        'L2Change_C',change,'estimatedFullStepTemporalError_C',2*change, ...
        'estimatedTemporalFractionOfFineL2',2*change/max(fine.metrics.L2_C,realmin), ...
        'note',['First-order backward-Euler estimate from dt versus dt/2 at fixed mesh. ' ...
                'It bounds observed temporal contamination, not an exact time-error theorem.']);
    if config.verbose
        fprintf('Finest dt/2: L2 %.6g C; time change %.3g C; estimated fraction %.3g.\n', ...
            fine.metrics.L2_C,change,timeRefinement.estimatedTemporalFractionOfFineL2);
    end
end
study=struct('benchmark','Warm saturated transient diffusion, mixed Dirichlet/Neumann', ...
    'version','v15_heat_diffusion_20260930','configuration',config, ...
    'caseResults',{caseResults},'comparison',comparison, ...
    'timeRefinement',timeRefinement,'totalWallSeconds',toc(started), ...
    'scope',['Analytic warm diffusion and hydrostatic-no-advection verification only; ' ...
             'no phase-front or gasket calibration is established by this patch.']);
if config.saveArtifacts
    folder=char(config.outputFolder); if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,'heat_diffusion_results.mat'),'study','-v7');
    writetable(comparison,fullfile(folder,'heat_diffusion_mesh_errors.csv'));
    if timeRefinement.performed
        timeTable=struct2table(timeRefinement.fine.metrics);
        timeTable.L2ChangeFromFullStep_C=timeRefinement.L2Change_C;
        timeTable.estimatedTemporalFraction=timeRefinement.estimatedTemporalFractionOfFineL2;
        writetable(timeTable,fullfile(folder,'heat_diffusion_time_check.csv'));
    end
    writeReport(study,fullfile(folder,'heat_diffusion_report.txt'));
end
end

function result=runCase(nx,dtSeconds,durationDays)
[x,z]=ndgrid(linspace(0,1,nx),linspace(0,1,nx)); xy=[x(:),z(:)];
node=reshape(1:nx^2,nx,nx); tri=zeros(2*(nx-1)^2,3); e=0;
for j=1:nx-1
    for i=1:nx-1
        a=node(i,j); b=node(i+1,j); c=node(i+1,j+1); d=node(i,j+1);
        % Positive orientation in mathematical (x,Z-down) coordinates.
        if mod(i+j,2)==0, pair=[a b d;b c d]; else, pair=[a b c;a c d]; end
        tri(e+(1:2),:)=pair; e=e+2;
    end
end
mesh=struct('pAll',xy,'tSoil',tri, ...
    'coordInfo',struct('sectionColumns','X [m], Z [m] positive down'));
opt=struct('mesh',mesh,'durationDays',durationDays,'outputEveryDays',durationDays, ...
    'dtInitial',dtSeconds,'dtMax',dtSeconds,'dtMin',dtSeconds/64, ...
    'topMin',2,'topMax',2,'topHTC',inf,'geothermalFlux',0, ...
    'hydroFar',true,'waterTableDepth',-5,'topSeepage',false,'tunnelDrain',false, ...
    'Ss',1e-18,'liningK',0,'advection',true,'enthalpyReference_Jkg',0, ...
    'soilViscosity',struct('enabled',false), ...
    'leakage',struct('jointFreeze',struct('enabled',false), ...
        'jointStorage',struct('enabled',false)), ...
    'mechanics',struct('enable',false),'plot',false,'verbose',false,'saveOutput',false, ...
    'solver',struct('tol',1e-11,'reuseAcrossSteps',true), ...
    'initial',struct('mode','fields','T',2+sin(pi*xy(:,2)/2), ...
        'h',xy(:,2)+5,'startDay',0));
started=tic; out=TH_2D_segment_v15(opt); wall=toc(started);
p=out.parameters;
capacity=(1-p.porosity)*p.rhoS*p.cs+p.rho*p.porosity*p.cw;
conductivity=(1-p.porosity)*p.ks+p.porosity*p.kw;
alpha=conductivity/capacity; timeSeconds=out.timeDay*p.day;
decay=exp(-alpha*(pi/2)^2*timeSeconds);
exact=2+decay*sin(pi*xy(:,2)/2);
L2=exactFieldError(out.T,mesh,decay);
top=abs(xy(:,2))<1e-12;
headError=max(abs(out.h-(xy(out.soilNodes,2)+5)));
initialExact=2+sin(pi*xy(:,2)/2);
% Ss is not silently removed by the solver. Bound its capacity correction.
relativeExtraCapacity=p.rho*p.cw*p.Ss*max(out.h)/capacity;
energyChange=max(1,abs(out.history(end,5)-out.history(1,5)));
metrics=struct('nx',nx,'nodes',size(xy,1),'elements',size(tri,1), ...
    'dx_m',1/(nx-1),'dtSeconds',dtSeconds,'finalDay',out.timeDay, ...
    'L2_C',L2,'nodalLinf_C',max(abs(out.T-exact)), ...
    'initialInterpolationL2_C',exactFieldError(initialExact,mesh,1), ...
    'alpha_m2_s',alpha,'exactDecayFactor',decay, ...
    'maxDarcyVelocity_m_s',max(abs(out.hydro.elemFlux(:))), ...
    'maxSubfaceFlow_m2_s',max(abs(out.hydro.subfaceFlux_m2s(:))), ...
    'maxHydrostaticHeadError_m',headError,'maxTopTemperatureError_C',max(abs(out.T(top)-2)), ...
    'maxIceVolumeFraction',max(out.state.thetaI), ...
    'relativeElasticHeatCapacity',relativeExtraCapacity, ...
    'cumulativeWaterDefect_m2',sum(out.history(:,8)), ...
    'cumulativeEnergyDefect_J_m',sum(out.history(:,9)), ...
    'relativeEnergyBudgetDefect',abs(sum(out.history(:,9)))/energyChange, ...
    'acceptedSteps',out.solverWork.acceptedSteps, ...
    'rejectedSteps',out.solverWork.rejectedSteps,'wallSeconds',wall);
result=struct('out',out,'exactNodalTemperature_C',exact,'metrics',metrics);
end

function value=exactFieldError(T,mesh,decay)
[N,w]=independentTriangleRule(); accumulated=0;
for e=1:size(mesh.tSoil,1)
    ids=mesh.tSoil(e,:); xy=mesh.pAll(ids,:);
    area=abs(det([xy(2,:)-xy(1,:);xy(3,:)-xy(1,:)]))/2;
    q=N*xy; difference=N*T(ids)-(2+decay*sin(pi*q(:,2)/2));
    accumulated=accumulated+area*sum(w.*difference.^2);
end
value=sqrt(accumulated);
end

function value=finiteElementDifference(T1,T2,mesh)
assert(numel(T1)==numel(T2)); [N,w]=independentTriangleRule(); accumulated=0;
for e=1:size(mesh.tSoil,1)
    ids=mesh.tSoil(e,:); xy=mesh.pAll(ids,:);
    area=abs(det([xy(2,:)-xy(1,:);xy(3,:)-xy(1,:)]))/2;
    difference=N*(T1(ids)-T2(ids));
    accumulated=accumulated+area*sum(w.*difference.^2);
end
value=sqrt(accumulated);
end

function [N,w]=independentTriangleRule()
% Positive Dunavant degree-5 rule; weights normalized to unit area.
a=.470142064105115; b=.059715871789770;
c=.101286507323456; d=.797426985353087;
N=[1/3 1/3 1/3;a a b;a b a;b a a;c c d;c d c;d c c];
w=[.225;repmat(.132394152788506,3,1);repmat(.125939180544827,3,1)];
end

function writeReport(study,path)
fid=fopen(path,'w'); assert(fid>=0,'Could not create benchmark report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',study.benchmark);
fprintf(fid,'Exact T = 2 + exp[-alpha*(pi/2)^2*t] sin(pi*z/2).\n');
fprintf(fid,'Unit square; top Dirichlet 2 C, sides/bottom insulated, h=z+5 m.\n');
fprintf(fid,'Warm saturated soil, no tunnel; advection enabled and viscosity disabled.\n');
fprintf(fid,'Mass and enthalpy use the v15 SCV control volumes.\n');
fprintf(fid,'L2 includes the finite-element interpolation error, using independent degree-5 integration.\n\n');
fprintf(fid,'nx     dx(m)       dt(s)       L2(C)        order       maxDarcy(m/s)\n');
for j=1:height(study.comparison)
    r=study.comparison(j,:);
    fprintf(fid,'%2d   %.8g   %.8g   %.8g   %.6g   %.4g\n', ...
        r.nx,r.dx_m,r.dtSeconds,r.L2_C,r.observedOrder,r.maxDarcyVelocity_m_s);
end
if study.timeRefinement.performed
    t=study.timeRefinement;
    fprintf(fid,'\nFinest mesh dt/2 L2 = %.9g C.\n',t.fine.metrics.L2_C);
    fprintf(fid,'L2(dt vs dt/2) = %.9g C.\n',t.L2Change_C);
    fprintf(fid,'Estimated full-step time error / fine L2 = %.6g.\n',t.estimatedTemporalFractionOfFineL2);
end
fprintf(fid,'\nTotal wall time %.4g s.\n',study.totalWallSeconds);
fprintf(fid,'Scope: %s\n',study.scope);
end
