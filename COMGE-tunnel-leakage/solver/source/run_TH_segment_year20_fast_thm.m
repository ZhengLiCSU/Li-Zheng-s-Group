function result=run_TH_segment_year20_fast_thm(user)
%RUN_TH_SEGMENT_YEAR20_FAST_THM Consecutive THM history with fast numerics.
% The physical model, nonlinear tolerance and admissible-state guards are
% those of run_TH_segment_year20_thm. This entry changes only the default
% mesh, positive SCV quadrature, time cap and Newton update strategy.
% March-April and September-October use at most two days; other months
% use at most four. CPU threads are scoped to this call and restored.
% Constitutive states and safeguarded Newton/LU tangents are reused; a
% failed stale-tangent line search immediately rebuilds the exact tangent.
% It retains all 20 years and the original installed mechanical reference.
%
% r=run_TH_segment_year20_fast_thm;
% r=run_TH_segment_year20_fast_thm(struct('resume',true));
% Options are the same as the strict annual driver. Nested model overrides
% can supply another physical case; its discretization accuracy must be
% checked separately. The supplied mesh protects the tunnel radius 7 m,
% topsoil depth 1 m and every original soil boundary and lining/joint.
% See efficiency_thm_20261001/README_efficiency_trial_zh.md for evidence
% and the distinction between an annual leakage error and monthly errors.

if nargin<1, user=struct(); end
assert(isstruct(user)&&isscalar(user),'THMFast:Config', ...
    'user must be a scalar struct.');
root=fileparts(mfilename('fullpath'));
model=struct('meshFile',fullfile(root,'example_mesh_thm_fast_R7_top1_d3.mat'), ...
    'dtMax',4*86400, ...
    'integration',struct('soilQuadratureOrder',2,'subdivisionLevel',0), ...
    'solver',struct('stepT',10,'stepH',250,'stepU',.01, ...
        'jacobianReuse',4,'jacobianRefreshRatio',.4,'reuseAcrossSteps',true, ...
        'reuseAcceptedResidual',true, ...
        'reuseConstitutiveState',true, ...
        'scvTangentAssembly','element', ...
        'dtMaxByMonth_days',[4 4 2 2 4 4 4 4 2 2 4 4]));
cfg=struct('maxYears',20,'resume',false, ...
    'outputFolder',fullfile(root,'segment_year20_thm_fast_v16_R7_top1_d3_illustrative'), ...
    'model',model,'numThreads',4);
for key=fieldnames(user)'
    assert(isfield(cfg,key{1}),'THMFast:Config','Unknown driver option: %s.',key{1});
    if strcmp(key{1},'model')
        assert(isstruct(user.model)&&isscalar(user.model),'THMFast:Config', ...
            'model must be a scalar struct.');
        cfg.model=combine(cfg.model,user.model);
    else
        cfg.(key{1})=user.(key{1});
    end
end
customMesh=isfield(user,'model')&& ...
    ((isfield(user.model,'mesh')&&~isempty(user.model.mesh)) ...
    ||isfield(user.model,'meshFile'));
if ~customMesh
    supplied=joint_supply_th_model(cfg.model);
    baseline=joint_supply_th_model();
    L=supplied.leakage; nj=numel(L.jointAnglesRad);
    conductance=L.jointConductance(:); widths=L.jointWidth_m(:);
    if isscalar(conductance), conductance=repmat(conductance,nj,1); end
    if isscalar(widths), widths=repmat(widths,nj,1); end
    assert(numel(conductance)==nj&&numel(widths)==nj,'THMFast:JointData', ...
        'Conductance and width must be scalar or one value per joint.');
    active=find(conductance>0);
    offset=L.jointAnglesRad(active)-baseline.leakage.jointAnglesRad(7);
    assert(all(abs(atan2(sin(offset),cos(offset)))<1e-10) ...
        &&all(abs(widths(active)-.02)<1e-12),'THMFast:CustomJointNeedsMesh', ...
        ['The supplied fast mesh resolves the default crown joint with 0.02 m ' ...
         'support. Supply a corresponding explicit mesh for other active joints.']);
end
threads=cfg.numThreads;
assert(isnumeric(threads)&&isreal(threads)&&isscalar(threads)&&isfinite(threads) ...
    &&threads>=1&&threads==fix(threads),'THMFast:Threads', ...
    'numThreads must be a positive integer.');
cfg=rmfield(cfg,'numThreads');
priorThreads=maxNumCompThreads;
cleanup=onCleanup(@()maxNumCompThreads(priorThreads)); %#ok<NASGU>
maxNumCompThreads(threads);
result=run_TH_segment_year20_thm(cfg);
result.executionThreads=threads;
execution=thm_fast_execution_metadata(result.lastYear,threads,priorThreads, ...
    result.elapsedWallTime_s);
save(fullfile(result.config.outputFolder,'fast_execution_metadata.mat'),'execution','-v7');
end

function a=combine(a,b)
for key=fieldnames(b)'
    name=key{1};
    if isfield(a,name)&&isstruct(a.(name))&&isstruct(b.(name))
        a.(name)=combine(a.(name),b.(name));
    else
        a.(name)=b.(name);
    end
end
end
