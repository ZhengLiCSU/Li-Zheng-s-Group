function study=run_thm_mesh_convergence(user)
%RUN_THM_MESH_CONVERGENCE Fresh histories on locally refined joint meshes.
% No displacement/installation/checkpoint transfer between meshes. The
% circular refiner changes local chords as well as element sizes, so this
% is joint-interface mesh/geometry convergence, NOT global soil-domain
% convergence. An input spatial field uses identical source coordinates and
% log values on every mesh, with no extrapolation of its physical support.
if nargin<1,user=struct();end
root=fileparts(mfilename('fullpath'));addpath(root);
c=struct('sourceMeshFile',fullfile(root,'example_mesh_repaired_interface.mat'), ...
    'outputFolder',fullfile(root,'leakage_research_20261001','mesh_histories'), ...
    'levels',[1 2 3],'maxYears',20,'model',struct(),'resume',false);
f=fieldnames(user);
for k=1:numel(f),assert(isfield(c,f{k}),'MeshStudy:Config','Unknown option.');c.(f{k})=user.(f{k});end
assert(isnumeric(c.levels)&&isvector(c.levels)&&numel(c.levels)>=2 ...
    &&all(c.levels>=1&c.levels<=5&c.levels==fix(c.levels)) ...
    &&all(diff(c.levels)>0),'MeshStudy:Levels','At least two increasing refiner levels in 1..5 are required.');
assert(~isfield(c.model,'mesh')&&~isfield(c.model,'meshFile'), ...
    'MeshStudy:MeshOverride','Supply sourceMeshFile, not a competing model.mesh/meshFile.');
if isfield(c.model,'initial')&&isfield(c.model.initial,'mode')
    assert(strcmp(c.model.initial.mode,'profile'),'MeshStudy:Initial', ...
        'Mesh histories require a shared profile, not state/installation transfer.');
end
S=load(c.sourceMeshFile);keyHash=fileHash(c.sourceMeshFile);
baseline=joint_contact_thm_model(combine(c.model,struct('mesh',S)));
L=baseline.leakage;active=find(L.jointConductance>0);
assert(~isempty(active),'MeshStudy:Joints','At least one existing defect is required.');
width=L.jointWidth_m;if isscalar(width),width=repmat(width,numel(L.jointAnglesRad),1);end
if ~isfolder(c.outputFolder),mkdir(c.outputFolder);end
metrics=cell(numel(c.levels),1);meshReports=cell(size(metrics));rows=struct([]);
for k=1:numel(c.levels)
    level=c.levels(k);folder=fullfile(c.outputFolder,sprintf('mesh_level_%d',level));
    if ~isfolder(folder),mkdir(folder);end
    meshFile=fullfile(folder,'local_mesh.mat');metadata=struct('sourceSHA256',keyHash, ...
        'level',level,'activeAnglesRad',L.jointAnglesRad(active),'jointWidths_m',width(active));
    if isfile(meshFile)
        d=load(meshFile,'mesh','meshReport','meshMetadata');
        assert(isequaln(d.meshMetadata,metadata),'MeshStudy:Existing', ...
            'Existing mesh belongs to another source/refinement; preserved.');
        mesh=d.mesh;meshReport=d.meshReport;
    else
        [mesh,meshReport]=refine_active_joint_mesh(S,L.jointAnglesRad(active),width(active),level);
        meshMetadata=metadata;save(meshFile,'mesh','meshReport','meshMetadata','-v7');
    end
    model=combine(c.model,struct('mesh',mesh));
    r=run_TH_segment_year20_thm(struct('outputFolder',fullfile(folder,'history'), ...
        'maxYears',c.maxYears,'model',model,'resume',c.resume));
    assert(r.completedYears==c.maxYears,'MeshStudy:Incomplete','Incomplete mesh history; evidence preserved.');
    metrics{k}=thm_seasonal_leakage_metrics(r.lastYear);meshReports{k}=meshReport;
    m=metrics{k};rows(k).level=level;rows(k).completedYears=r.completedYears;
    rows(k).nodes=size(mesh.pAll,1);rows(k).soilTriangles=size(mesh.tSoil,1);
    rows(k).annual_L_m=m.annual_L_m;rows(k).peakDaily_L_m=m.peakDaily_L_m;
    rows(k).firstActiveDay=m.firstActiveDay;rows(k).lastActiveDay=m.lastActiveDay;
    save(fullfile(c.outputFolder,'mesh_study_progress.mat'),'c','metrics','meshReports','rows','-v7');
end
comparison=cellfun(@(m)compare_thm_seasonal_metrics(m,metrics{end}),metrics,'UniformOutput',false);
study=struct('complete',true,'isFull20HistoryStudy',c.maxYears==20,'config',c, ...
    'summary',struct2table(rows),'metrics',{metrics},'meshReports',{meshReports}, ...
    'comparisonToLastMesh',{comparison},'scope', ...
    'Independent local joint mesh/geometry histories with fixed physics; not global soil convergence or experimental validation.');
save(fullfile(c.outputFolder,'mesh_convergence_report.mat'),'study','-v7');
writetable(study.summary,fullfile(c.outputFolder,'mesh_convergence_summary.csv'));
end
function a=combine(a,b)
f=fieldnames(b);
for k=1:numel(f)
    n=f{k};if isfield(a,n)&&isstruct(a.(n))&&isstruct(b.(n)),a.(n)=combine(a.(n),b.(n));else,a.(n)=b.(n);end
end
end
function digest=fileHash(file)
fid=fopen(file,'rb');assert(fid>0,'MeshStudy:Source','Cannot open source mesh.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes(:),'int8'));raw=typecast(md.digest(),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
end
