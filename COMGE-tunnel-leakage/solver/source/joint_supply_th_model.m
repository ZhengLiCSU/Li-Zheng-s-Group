function model=joint_supply_th_model(overrides)
%JOINT_SUPPLY_TH_MODEL Verified-physics input builder for the v14 entry.
% Soil-water viscosity, equilibrium channel storage and joint icing are
% enabled by default. Nested overrides preserve unspecified baseline data.
% Explicit nonempty mesh or meshFile inputs are returned unchanged.
%
% Without an explicit mesh, the supplied source mesh is locally refined at
% the default active crown joint (Joint 7, 0.02 m support, Level 1). The
% derived mesh is cached beside this function with an absolute file path.
% It includes source/refiner SHA-256 fingerprints to reject stale caches.
% A custom active joint location or support width requires a corresponding
% explicit mesh; the default crown refinement is not claimed to resolve it.
% The v12 baseline and its source mesh are not modified.
if nargin<1, overrides=struct(); end
assert(isstruct(overrides)&&isscalar(overrides), ...
    'joint_supply_th_model:InvalidOverrides','Overrides must be a scalar struct.');
root=fileparts(mfilename('fullpath'));
physical=struct('soilViscosity',struct('enabled',true,'referenceTemperature_C',20), ...
    'leakage',struct('jointStorage',struct('enabled',true,'tortuosity',1)));
model=joint_freeze_th_model(combine(physical,overrides));
if isfield(overrides,'mesh')&&~isempty(overrides.mesh)
    return
end
if isfield(overrides,'meshFile')
    assert((ischar(overrides.meshFile)|| ...
        (isstring(overrides.meshFile)&&isscalar(overrides.meshFile))) ...
        &&~isempty(char(overrides.meshFile)), ...
        'joint_supply_th_model:InvalidMeshFile', ...
        'An explicit meshFile must be a nonempty character vector or scalar string.');
    return
end

baseline=joint_freeze_th_model();
targetAngle=baseline.leakage.jointAnglesRad(7); targetWidth=.02; level=1;
L=model.leakage; angles=double(L.jointAnglesRad(:)); nj=numel(angles);
c=jointVector(L.jointConductance,nj,'jointConductance');
widths=jointVector(L.jointWidth_m,nj,'jointWidth_m');
active=find(c>0);
distance=abs(atan2(sin(angles(active)-targetAngle),cos(angles(active)-targetAngle)));
assert(all(distance<1e-10)&&all(abs(widths(active)-targetWidth)<1e-12), ...
    'joint_supply_th_model:CustomJointNeedsMesh', ...
    ['Automatic refinement resolves only the default crown joint at %.9g rad ' ...
     'with 0.02 m support. Supply mesh or meshFile for custom active joints or widths.'], ...
    targetAngle);
source=fullfile(root,'example_mesh_repaired_interface.mat');
refiner=fullfile(root,'refine_active_joint_mesh.m');
assert(isfile(source)&&isfile(refiner),'joint_supply_th_model:MissingMeshSource', ...
    'The default source mesh and refinement function must exist beside this file.');
meshFile=fullfile(root,'example_mesh_joint7_refined_level1.mat');
signature=struct('schemaVersion',1,'sourceSHA256',fileHash(source), ...
    'refinerSHA256',fileHash(refiner),'activeAngleRad',targetAngle, ...
    'jointWidth_m',targetWidth,'level',level);
validCache=false;
if isfile(meshFile)
    variables=whos('-file',meshFile);
    if all(ismember({'jointRefinementCacheSignature','jointRefinementReport'}, ...
            {variables.name}))
        cached=load(meshFile,'jointRefinementCacheSignature','jointRefinementReport');
        validCache=isequaln(cached.jointRefinementCacheSignature,signature);
        if validCache
            r=cached.jointRefinementReport;
            validCache=r.level==level&&r.negativeInterfaceHydraulicEdges==0 ...
                &&r.retainedUnusedNodes==0 ...
                &&all(r.actualBandEdgeLengthRange_m(:,2)<=targetWidth/2^level*(1+1e-6));
        end
    end
end
if ~validCache
    S=load(source);
    [Sfine,report]=refine_active_joint_mesh(S,targetAngle,targetWidth,level);
    assert(report.negativeInterfaceHydraulicEdges==0&&report.retainedUnusedNodes==0, ...
        'joint_supply_th_model:InvalidRefinement','Generated mesh failed interface checks.');
    Sfine.jointRefinementCacheSignature=signature;
    temporary=[tempname(root) '.mat'];
    clean=onCleanup(@()removeTemporary(temporary)); %#ok<NASGU>
    save(temporary,'-struct','Sfine','-v7');
    [success,message]=movefile(temporary,meshFile,'f');
    assert(success,'joint_supply_th_model:CacheWriteFailed', ...
        'Cannot store the derived mesh: %s',message);
end
model.meshFile=meshFile;
end

function a=combine(a,b)
names=fieldnames(b);
for k=1:numel(names)
    name=names{k};
    if isfield(a,name)&&isstruct(a.(name))&&isstruct(b.(name)) ...
            &&isscalar(a.(name))&&isscalar(b.(name))
        a.(name)=combine(a.(name),b.(name));
    else
        a.(name)=b.(name);
    end
end
end

function x=jointVector(x,n,name)
assert(isnumeric(x)&&isreal(x)&&(isvector(x)||isempty(x)) ...
    &&all(isfinite(x(:)))&&(isscalar(x)||numel(x)==n), ...
    'joint_supply_th_model:InvalidJointParameter', ...
    '%s must be finite, scalar or one value per joint.',name);
x=double(x(:)); if isscalar(x),x=repmat(x,n,1); end
end

function digest=fileHash(file)
fid=fopen(file,'rb');
assert(fid>0,'joint_supply_th_model:CannotReadSource','Cannot read source: %s',file);
clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,inf,'*uint8'); md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes(:),'int8')); raw=typecast(md.digest(),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
end

function removeTemporary(file)
if isfile(file), delete(file); end
end
