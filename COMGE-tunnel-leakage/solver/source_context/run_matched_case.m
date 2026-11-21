function run_matched_case(root,caseid)
folder=fullfile(root,'cases',caseid);
assert(isfolder(folder)&&~isfile(fullfile(folder,'started.json')),'RFMesh:Once','Case must be fresh');
addpath(fullfile(root,'source'));maxNumCompThreads(4);
t=tic;while ~isfile(fullfile(folder,'observer_ready.json')),assert(toc(t)<60,'RFMesh:Observer','Observer did not acquire process handle');pause(.1);end
plan=jsondecode(fileread(fullfile(root,'plan.json')));cfg=plan.cases(strcmp({plan.cases.id},caseid));assert(isscalar(cfg));
old=load(fullfile(root,'baseline_year1.mat'),'out');p=old.out.parameters;
input=load(fullfile(root,'inputs',[caseid '.mat']));
p.meshFile=fullfile(root,'meshes',[cfg.mesh '.mat']);p.durationDays=1825;p.initial.startDay=0;
p.soilKsatField=struct('coordinates_m',input.soil_coordinates_m,'logMultiplier',input.soil_log_multiplier,'quadratureCoordinates_m',input.quadrature_coordinates_m,'quadratureLogMultiplier',input.quadrature_log_multiplier);
p.leakage.jointFreeze.quadraturePoints=64;p.leakage.jointFreeze.flowWeightedRadiusProfile_m=input.joint_radius_profile_m;
p.outputPrefix=fullfile(folder,'full1825');p.verbose=false;p.plot=false;p.saveOutput=true;
assert(p.thm.enabled&&~p.mechanics.enable&&strcmp(p.initial.mode,'profile'));assert(isequal(size(input.joint_radius_profile_m),[8,64]));
save(fullfile(folder,'executed_parameters.mat'),'p','cfg','-v7');
writejson(fullfile(folder,'started.json'),struct('utc',char(datetime('now','TimeZone','UTC')),'caseid',caseid,'days',1825));
wall=tic;cpu=cputime;
try
 out=TH_2D_segment_v16(p);
 assert(abs(out.timeDay-1825)<1e-8&&all(isfinite(out.history(:))));assert(max(out.history(:,11))<=p.solver.tol*1.001);
 writejson(fullfile(folder,'completed.json'),struct('solver_returned',true,'days',out.timeDay,'wall_s',toc(wall),'script_cpu_s',cputime-cpu,'steps',size(out.history,1),'max_scaled_residual',max(out.history(:,11))));
 fprintf('MATCHED_RF_MESH_COMPLETED %s %.12g days\n',caseid,out.timeDay);
catch ME
 writejson(fullfile(folder,'failure.json'),struct('identifier',ME.identifier,'message',ME.message,'stack',ME.stack,'wall_s',toc(wall),'script_cpu_s',cputime-cpu));rethrow(ME)
end
end
function writejson(p,x)
assert(~isfile(p));f=fopen(p,'w','n','UTF-8');assert(f>=0);cleanup=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(x,PrettyPrint=true));
end
