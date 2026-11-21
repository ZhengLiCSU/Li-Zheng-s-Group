function report=run_joint_supply_benchmark(outputDir)
%RUN_JOINT_SUPPLY_BENCHMARK Warm saturated linear soil/joint resistance.
% This is an independent steady Laplace hydraulic benchmark, not a frozen
% TH solution or a year-20 prediction. It holds K uniform and saturated,
% uses nominal 20 C reference conductances, and excludes storage, rainfall,
% soil/joint phase change, viscosity changes and thermal advection.
% Side and bottom total head follow the 5 m water table. The surface and
% concrete body are no-flow. Only joint 7 of the eight-segment ring leaks.
% Original and two locally refined meshes use consistent finite-band Robin
% integration. No negative FEM weights are clipped.
%
% REPORT=RUN_JOINT_SUPPLY_BENCHMARK(OUTPUTDIR) writes tables and MAT details.
% Each mesh/K/C case is a sparse linear solve. Refined meshes are generated
% from the original mesh, with fresh independent steady fields.
root=fileparts(mfilename('fullpath'));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(root,'joint_supply_study_20260930','warm_linear_benchmark');
end
outputDir=char(outputDir);
assert(~isempty(outputDir),'Output directory cannot be empty.');
if ~isfolder(outputDir), mkdir(outputDir); end
meshFile=fullfile(root,'example_mesh_repaired_interface.mat');
original=load(meshFile);
base=joint_only_th_model(); L=base.leakage;
assert(numel(L.jointAnglesRad)==8,'This benchmark requires the eight-joint example.');
activeJoint=7; width=L.jointWidth_m;
if ~isscalar(width), width=width(activeJoint); end
levels=[0 1 2]; permeabilities=[1e-7 1e-9 1e-10];
conductances=[2e-11 2e-10 2e-9];
waterTableDepth=5; secondsPerYear=365*86400;
segmentBounds=sort(mod(L.jointAnglesRad(:),2*pi));
rows=struct([]); details=struct([]); meshes=cell(numel(levels),1);
meshReports=cell(numel(levels),1); record=0;
for meshIndex=1:numel(levels)
    level=levels(meshIndex);
    if level==0
        S=original;
        meshReport=struct('level',0,'method','original repaired interface mesh', ...
            'refinedNodes',size(S.pAll,1),'refinedSoilTriangles',size(S.tSoil,1), ...
            'refinedInterfaceEdges',size(S.soilTunnelCoupleTable,1));
    else
        [S,meshReport]=refine_active_joint_mesh(original, ...
            L.jointAnglesRad(activeJoint),width,level);
    end
    meshes{meshIndex}=S; meshReports{meshIndex}=meshReport;
    nodes=unique(S.tSoil(:)); ns=numel(nodes); xy=S.pAll(nodes,:);
    map=zeros(size(S.pAll,1),1); map(nodes)=1:ns;
    local=map(S.tSoil);
    [A0,~,geometricWeight]=soilStiffness(xy,local);
    bounds=min(xy,[],1); maxima=max(xy,[],1); tol=1e-8;
    far=find(abs(xy(:,1)-bounds(1))<tol ...
        |abs(xy(:,1)-maxima(1))<tol|abs(xy(:,2)-maxima(2))<tol);
    free=setdiff((1:ns)',far);
    assert(~isempty(far)&&~isempty(free),'No far-head boundary or interior soil nodes.');
    Q=build_joint_drain_quadrature(S,nodes,L,segmentBounds);
    take=Q.jointID==activeJoint;
    phi=Q.phi(take,:); weights=Q.weight(take);
    assert(abs(sum(weights)-1)<1e-12,'Active joint weights are not normalized.');
    depth=xy(:,2)-min(S.pAll(:,2));
    interfaceNodes=unique(S.soilTunnelCoupleTable(:,2:3));
    innerNodes=S.tunnelInnerRingGlobal(:);
    innerRadius=mean(hypot(S.pAll(innerNodes,1)-S.xc,S.pAll(innerNodes,2)-S.yc));
    outerRadius=hypot(S.pAll(interfaceNodes,1)-S.xc,S.pAll(interfaceNodes,2)-S.yc);
    thickness=outerRadius-innerRadius;
    drop=zeros(ns,1);
    drop(map(interfaceNodes))=-thickness.* ...
        (S.pAll(interfaceNodes,2)-S.yc)./outerRadius;
    interiorHead=L.interiorHead;
    entry=L.jointEntryHead;
    if isempty(entry), entry=0; elseif ~isscalar(entry), entry=entry(activeJoint); end
    referenceNodal=depth-waterTableDepth-interiorHead+drop;
    referenceQuad=phi*referenceNodal-entry;
    assert(all(referenceQuad>0),'The example drain must lie below its reference head.');
    H0=weights'*referenceQuad;
    W=spdiags(weights,0,numel(weights),numel(weights));
    robin=phi'*W*phi;
    unitLoad=phi'*weights;
    A0ff=A0(free,free);
    [~,cholFlag]=chol(A0ff);
    assert(cholFlag==0,'Free soil Laplace matrix is not positive definite.');
    shapeResistance=unitLoad(free)'*(A0ff\unitLoad(free));
    assert(isfinite(shapeResistance)&&shapeResistance>0,'Invalid soil spreading resistance.');
    supportedNodes=nodes(find(full(sum(abs(phi),1))>0)); %#ok<FNDSB>
    interfacePairs=S.soilTunnelCoupleTable(:,2:3);
    touched=any(ismember(interfacePairs,supportedNodes),2);
    lengths=hypot(S.pAll(interfacePairs(:,1),1)-S.pAll(interfacePairs(:,2),1), ...
        S.pAll(interfacePairs(:,1),2)-S.pAll(interfacePairs(:,2),2));
    maximumBandEdge=max(lengths(touched));
    fprintf('Warm linear benchmark: mesh level %d, %d soil nodes, support %.6g m.\n', ...
        level,ns,Q.supportLength(activeJoint));
    for K=permeabilities
        A=K*A0; Gsoil=K/shapeResistance;
        for C=conductances
            % Drawdown from the undisturbed hydrostatic total head avoids
            % cancellation of a constant potential in the FE operator.
            matrix=A+C*robin;
            rhs=C*phi'*(weights.*referenceQuad);
            drawdown=zeros(ns,1);
            drawdown(free)=matrix(free,free)\rhs(free);
            head=depth-waterTableDepth-drawdown;
            actualQuad=referenceQuad-phi*drawdown;
            assert(all(actualQuad>=-1e-9), ...
                ['Chosen case activates reverse Robin flow. A unilateral ' ...
                 'active-set solve is required beyond this linear benchmark.']);
            rateQuad=C*weights.*actualQuad;
            qNode=phi'*rateQuad; rate=sum(rateQuad);
            hydraulicResidual=-A*drawdown+qNode;
            farReaction=zeros(ns,1); farReaction(far)=-hydraulicResidual(far);
            waterDefect=sum(qNode+farReaction);
            hydroRate=C*H0;
            seriesRate=H0/(1/C+1/Gsoil);
            actualHead=weights'*actualQuad;
            record=record+1;
            caseRow=struct('meshLevel',level,'soilNodes',ns, ...
                'soilTriangles',size(S.tSoil,1), ...
                'jointSupportLength_m',Q.supportLength(activeJoint), ...
                'maximumSupportedInterfaceEdge_m',maximumBandEdge, ...
                'K_m_per_s',K,'C_m_per_s',C,'referenceDrivingHead_m',H0, ...
                'actualDrivingHead_m',actualHead,'Q_m2_per_s',rate, ...
                'Q_L_per_day_m',rate*86400*1000, ...
                'steadyEquivalentAnnual_L_per_m',rate*secondsPerYear*1000, ...
                'C_times_staticHead_m2_per_s',hydroRate, ...
                'staticHeadRelativeDifference_percent',100*abs(rate-hydroRate)/hydroRate, ...
                'uniformBandSoilConductance_m_per_s',Gsoil,'Lambda',C/Gsoil, ...
                'uniformBandSeriesRate_m2_per_s',seriesRate, ...
                'seriesRateRelativeDifference_percent',100*abs(rate-seriesRate)/seriesRate, ...
                'minimumTransmissivity_m_per_s',min(K*geometricWeight), ...
                'negativeTransmissivityEdges',nnz(geometricWeight< -1e-12), ...
                'massDefect_m2_per_s',waterDefect, ...
                'maximumFreeResidual_m2_per_s',norm(hydraulicResidual(free),inf), ...
                'minimumQuadratureDrivingHead_m',min(actualQuad));
            detailRow=struct('meshIndex',meshIndex,'meshLevel',level, ...
                'K_m_per_s',K,'C_m_per_s',C,'soilNodes',nodes, ...
                'head_m',head,'drawdown_m',drawdown,'nodalDrain_m2_per_s',qNode, ...
                'farReaction_m2_per_s',farReaction,'quadrature',Q, ...
                'activeJoint',activeJoint,'actualQuadratureDrivingHead_m',actualQuad, ...
                'referenceQuadratureDrivingHead_m',referenceQuad);
            if record==1
                rows=caseRow; details=detailRow;
            else
                rows(record)=caseRow; details(record)=detailRow; %#ok<AGROW>
            end
        end
    end
end
results=struct2table(rows);
convergence=convergenceTable(results);
report=struct('kind','warm saturated steady linear hydraulic benchmark', ...
    'meshFile',meshFile,'outputDir',outputDir,'waterTableDepth_m',waterTableDepth, ...
    'activeJoint',activeJoint,'results',results,'convergence',convergence, ...
    'meshReports',{meshReports}, ...
    'soilConductanceDefinition', ...
    'G=1/(w''*A_ff^(-1)*w), w is a normalized uniform-band nodal load', ...
    'interpretation', ...
    ['These steady warm rates are not frozen TH or year-20 results. ' ...
     'The scalar series formula uses a prescribed uniform load shape; ' ...
     'finite-band Robin redistribution can differ from it. Refinement ' ...
     'slightly changes the polygonal tunnel geometry.']);
writetable(results,fullfile(outputDir,'warm_linear_supply.csv'));
writetable(convergence,fullfile(outputDir,'warm_linear_mesh_convergence.csv'));
save(fullfile(outputDir,'warm_linear_supply.mat'),'report','details','meshes','-v7.3');
disp(results(:,{'meshLevel','K_m_per_s','C_m_per_s','Q_m2_per_s', ...
    'staticHeadRelativeDifference_percent','Lambda','actualDrivingHead_m'}));
end

function [A,edges,weight]=soilStiffness(xy,tri)
% Standard exact T3 Laplace stiffness; no permeability or Robin terms.
a=xy(tri(:,1),:); b=xy(tri(:,2),:); c=xy(tri(:,3),:);
twice=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2)) ...
    -(c(:,1)-a(:,1)).*(b(:,2)-a(:,2));
assert(all(abs(twice)>1e-12),'Degenerate soil triangle in the benchmark.');
area=abs(twice)/2;
gradient=cat(3,[b(:,2)-c(:,2),c(:,1)-b(:,1)]./twice, ...
    [c(:,2)-a(:,2),a(:,1)-c(:,1)]./twice, ...
    [a(:,2)-b(:,2),b(:,1)-a(:,1)]./twice);
pairs=[tri(:,[1 2]);tri(:,[2 3]);tri(:,[3 1])];
value=-[area.*sum(gradient(:,:,1).*gradient(:,:,2),2); ...
    area.*sum(gradient(:,:,2).*gradient(:,:,3),2); ...
    area.*sum(gradient(:,:,3).*gradient(:,:,1),2)];
[edges,~,edgeMap]=unique(sort(pairs,2),'rows');
weight=accumarray(edgeMap,value,[size(edges,1) 1]);
nf=size(edges,1); ns=size(xy,1);
inc=sparse([edges(:,1);edges(:,2)],[(1:nf)';(1:nf)'], ...
    [ones(nf,1);-ones(nf,1)],ns,nf);
A=inc*spdiags(weight,0,nf,nf)*inc';
end

function T=convergenceTable(results)
% Compare like K/C pairs, without treating the finest mesh as exact truth.
pairs=unique([results.K_m_per_s results.C_m_per_s],'rows');
rows=struct([]); record=0;
for pair=1:size(pairs,1)
    mask=results.K_m_per_s==pairs(pair,1)&results.C_m_per_s==pairs(pair,2);
    R=sortrows(results(mask,:),'meshLevel');
    for k=2:height(R)
        record=record+1;
        caseRow=struct('K_m_per_s',pairs(pair,1),'C_m_per_s',pairs(pair,2), ...
            'coarseLevel',R.meshLevel(k-1),'fineLevel',R.meshLevel(k), ...
            'coarseRate_m2_per_s',R.Q_m2_per_s(k-1),'fineRate_m2_per_s',R.Q_m2_per_s(k), ...
            'rateChangeRelativeToFine_percent', ...
            100*abs(R.Q_m2_per_s(k)-R.Q_m2_per_s(k-1))/abs(R.Q_m2_per_s(k)), ...
            'soilConductanceChangeRelativeToFine_percent', ...
            100*abs(R.uniformBandSoilConductance_m_per_s(k) ...
            -R.uniformBandSoilConductance_m_per_s(k-1)) ...
            /abs(R.uniformBandSoilConductance_m_per_s(k)));
        if record==1
            rows=caseRow;
        else
            rows(record)=caseRow; %#ok<AGROW>
        end
    end
end
T=struct2table(rows);
end
