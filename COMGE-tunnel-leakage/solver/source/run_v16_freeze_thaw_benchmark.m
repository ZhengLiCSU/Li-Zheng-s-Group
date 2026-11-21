function study=run_v16_freeze_thaw_benchmark(options)
%RUN_V16_FREEZE_THAW_BENCHMARK Actual v16 against independent phase fronts.
% Both cooling and warming are solved with the production residual/Jacobian
% on a 2-D strip. Isolated TH verification deliberately disables mechanics
% and joints; the existing coupled/contact tests verify those components.
% Physical benchmark constants differ from the illustrative tunnel case.
if nargin<1,options=struct();end
root=fileparts(mfilename('fullpath'));
c=struct('nz',[65 129 257],'dtSeconds',300,'durationDays',1,'startDay',1, ...
    'depth_m',3,'width_m',.02,'checkTimeRefinement',true, ...
    'saveArtifacts',true,'verbose',true, ...
    'outputFolder',fullfile(root,'physics_validation_20261001','phase_front'));
for name=fieldnames(options)'
    assert(isfield(c,name{1}),'PhaseBenchmark:Option','Unknown option: %s',name{1});
    c.(name{1})=options.(name{1});
end
assert(isvector(c.nz)&&all(c.nz>=3)&&all(c.nz==fix(c.nz))&&all(diff(c.nz)>0));
assert(c.dtSeconds>0&&c.durationDays>0&&c.startDay>0&&c.depth_m>0&&c.width_m>0);
cases=cell(2,numel(c.nz));references=cell(2,1);records=struct([]);n=0;
timeRefinement=cell(2,1);timer=tic;
for mode=1:2
    if mode==1,Ts=-4;Ti=2;label='freezing';else,Ts=4;Ti=-2;label='thawing';end
    ref=freeze_thaw_similarity_reference(Ts,Ti);references{mode}=ref;
    for j=1:numel(c.nz)
        cases{mode,j}=runCase(c.nz(j),c.dtSeconds,c,ref,label);
        n=n+1;
        if n==1,records=cases{mode,j}.metrics;
        else,records(n,1)=cases{mode,j}.metrics;end %#ok<AGROW>
        if c.verbose,disp(records(n));end
    end
    timeRefinement{mode}=struct('performed',false);
    if c.checkTimeRefinement
        fine=runCase(c.nz(end),c.dtSeconds/2,c,ref,label);
        previous=cases{mode,end};
        timeRefinement{mode}=struct('performed',true,'fine',fine, ...
            'maxTemperatureChange_C',max(abs(fine.out.T-previous.out.T)), ...
            'frontChange_m',abs(fine.metrics.computedFront_m-previous.metrics.computedFront_m));
    end
end
comparison=struct2table(records);
study=struct('configuration',c,'references',{references},'cases',{cases}, ...
    'comparison',comparison,'timeRefinement',{timeRefinement},'elapsed_s',toc(timer), ...
    'scope',['Independent moving freezing/thawing fronts for the saturated, constant-k TH limit. ' ...
    'Does not verify frost damage, ice lenses, or field/gasket parameters.']);
if c.saveArtifacts
    if ~isfolder(c.outputFolder),mkdir(c.outputFolder);end
    save(fullfile(c.outputFolder,'freeze_thaw_benchmark.mat'),'study','-v7');
    writetable(comparison,fullfile(c.outputFolder,'phase_front_errors.csv'));
end
end

function result=runCase(nz,dt,c,reference,label)
% Three columns give internal hydraulic DOFs as well as prescribed far heads.
[x,z]=ndgrid(linspace(0,c.width_m,3),linspace(0,c.depth_m,nz));xy=[x(:),z(:)];
node=reshape(1:numel(x),3,nz);tri=zeros(4*(nz-1),3);e=0;
for j=1:nz-1
    for i=1:2
        a=node(i,j);b=node(i+1,j);cc=node(i+1,j+1);d=node(i,j+1);
        if mod(i+j,2)==0,pair=[a b d;b cc d];else,pair=[a b cc;a cc d];end
        tri(e+(1:2),:)=pair;e=e+2;
    end
end
mesh=struct('pAll',xy,'tSoil',tri);
p=reference.parameters;
initial=deval(reference.solution,xy(:,2)'/(2*sqrt(reference.alpha_m2_s*c.startDay*86400)),1)';
opt=struct('mesh',mesh,'durationDays',c.durationDays,'dtInitial',dt,'dtMax',dt, ...
    'dtMin',dt/64,'outputEveryDays',c.durationDays/4, ...
    'porosity',p.porosity,'thetaR',p.thetaR,'vgAlpha',p.vgAlpha,'vgN',p.vgN, ...
    'kw',2,'ki',2,'ks',2,'ka',2,'Ksat',1e-7,'Ss',1e-10, ...
    'topMin',reference.surfaceTemperature_C,'topMax',reference.surfaceTemperature_C, ...
    'topHTC',inf,'geothermalFlux',0,'waterTableDepth',-5,'hydroFar',true, ...
    'topSeepage',false,'tunnelDrain',false,'liningK',0,'advection',true, ...
    'soilViscosity',struct('enabled',false), ...
    'thm',struct('enabled',false),'mechanics',struct('enable',false), ...
    'leakage',struct('jointFreeze',struct('enabled',false), ...
        'jointStorage',struct('enabled',false)), ...
    'initial',struct('mode','fields','T',initial,'h',xy(:,2)+5,'startDay',c.startDay), ...
    'solver',struct('tol',1e-11,'maxIter',60,'stepT',2,'reuseAcrossSteps',true), ...
    'plot',false,'verbose',false,'saveOutput',false);
timer=tic;out=TH_2D_segment_v16(opt);elapsed=toc(timer);
endTime=out.timeDay*86400;
scale=2*sqrt(reference.alpha_m2_s*endTime);
exact=deval(reference.solution,xy(:,2)'/scale,1)';
Tmean=mean(reshape(out.T,3,nz),1)';zz=z(1,:)';
id=find(Tmean(1:end-1).*Tmean(2:end)<=0,1);
assert(~isempty(id),'PhaseBenchmark:Front','The computed phase front is missing.');
front=zz(id)-Tmean(id)*(zz(id+1)-zz(id))/(Tmean(id+1)-Tmean(id));
exactFront=reference.zeroTemperatureEta*scale;
L2=independentL2(out.T,mesh,reference,scale,c.width_m);
headError=max(abs(out.h-(xy(out.soilNodes,2)+5)));
initialEnergy=(1-p.porosity)*p.rhoS*p.cs+p.rho*p.porosity*p.cw;
relativeElasticCapacity=out.parameters.rho*out.parameters.cw*out.parameters.Ss ...
    *max(out.h)/initialEnergy;
metrics=struct('mode',label,'nz',nz,'dz_m',c.depth_m/(nz-1),'dtSeconds',dt, ...
    'finalDay',out.timeDay,'temperatureL2_C',L2,'temperatureNodalLinf_C',max(abs(out.T-exact)), ...
    'computedFront_m',front,'referenceFront_m',exactFront, ...
    'frontError_m',abs(front-exactFront),'frontRelativeError',abs(front-exactFront)/exactFront, ...
    'relativeElasticCapacity',relativeElasticCapacity,'maxHeadError_m',headError, ...
    'maxDarcyVelocity_m_s',max(abs(out.hydro.elemFlux(:))), ...
    'maxWaterDefect_m2',max(abs(out.history(:,8))), ...
    'maxEnergyDefect_J_m',max(abs(out.history(:,9))), ...
    'cumulativeEnergyDefect_J_m',sum(out.history(:,9)), ...
    'energyChange_J_m',out.history(end,5)-out.history(1,5), ...
    'acceptedSteps',out.solverWork.acceptedSteps,'rejectedSteps',out.solverWork.rejectedSteps, ...
    'wallSeconds',elapsed);
result=struct('out',out,'referenceNodalTemperature_C',exact,'metrics',metrics);
end

function L2=independentL2(T,mesh,reference,scale,width)
a=.470142064105115;b=.059715871789770;c=.101286507323456;d=.797426985353087;
N=[1/3 1/3 1/3;a a b;a b a;b a a;c c d;c d c;d c c];
w=[.225;repmat(.132394152788506,3,1);repmat(.125939180544827,3,1)];total=0;
for j=1:size(mesh.tSoil,1)
    ids=mesh.tSoil(j,:);xy=mesh.pAll(ids,:);q=N*xy;
    exact=deval(reference.solution,q(:,2)'/scale,1)';
    area=abs(det([xy(2,:)-xy(1,:);xy(3,:)-xy(1,:)]))/2;
    total=total+area*sum(w.*(N*T(ids)-exact).^2);
end
L2=sqrt(total/(width*max(mesh.pAll(:,2))));
end
