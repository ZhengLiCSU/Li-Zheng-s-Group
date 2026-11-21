function [reduced,pReduced] = reduce_premelting_cyclic_symmetry(fullModel,p)
%REDUCE_PREMELTING_CYCLIC_SYMMETRY Exact eightfold projection for NOMINAL runs.
% Uniform baths/ports and identical aperture profiles are required. General
% random samples use the full model. All operators retain one-sector energy.
arguments
    fullModel (1,1) struct
    p (1,1) struct
end
M=fullModel.mesh; ids=find(M.nodeSegmentID==1); nt=numel(ids); ns=numel(fullModel.capacity); nq=fullModel.nq;
Pt=sparse(ns,nt); Pu=sparse(2*ns,2*nt);
for j=1:8
    group=find(M.nodeSegmentID==j);
    assert(numel(group)==nt,'PhaseSwitch:Symmetry','Segment node counts differ.');
    angle=(j-1)*pi/4; R=[cos(angle),-sin(angle);sin(angle),cos(angle)];
    expected=(M.p(ids,:)-M.center)*R'+M.center;
    assert(max(abs(M.p(group,:)-expected),[],'all')<1e-10,'PhaseSwitch:Symmetry','Geometry lacks cyclic symmetry.');
    Pt(group,:)=speye(nt);
    for k=1:nt, Pu(2*group(k)+[-1,0],2*k+[-1,0])=R; end
end
cellRows=1:nq; faceRows=1:nq+1; jointMap=repmat(speye(nq),8,1);
B=fullModel.gapOperator*Pu;
assert(norm(B-jointMap*B(cellRows,:),inf)<1e-10,'PhaseSwitch:Symmetry','Joint operator is not cyclic.');
assert(norm(Pu'*(fullModel.gapOperator'*diag(fullModel.weight))*jointMap/8-B(cellRows,:)'*diag(fullModel.weight(cellRows)),inf)<1e-10, ...
    'PhaseSwitch:Symmetry','Mechanical port projection failed.');
reduced=fullModel; reduced.weight=fullModel.weight(cellRows);
reduced.capacity=Pt'*fullModel.capacity/8;
reduced.stiffness=Pu'*fullModel.stiffness*Pu/8;
reduced.thermalStiffness=Pt'*fullModel.thermalStiffness*Pt/8;
reduced.thermalLoad=Pu'*fullModel.thermalLoad*Pt/8;
reduced.gapOperator=B(cellRows,:);
reduced.heatExchange=fullModel.heatExchange(cellRows,:)*Pt;
reduced.exchangeRowSum=sum(reduced.heatExchange,2);
reduced.exchangeColumnSum=sum(reduced.heatExchange,1)';
names={'bathConductance','innerBathConductance','outerBathConductance'};
for k=1:numel(names), reduced.(names{k})=Pt'*fullModel.(names{k})/8; end
reduced.faces=fullModel.faces(faceRows,:); reduced.faceDistance=fullModel.faceDistance(faceRows); reduced.boundaryID=fullModel.boundaryID(faceRows);
reduced.divergence=fullModel.divergence(cellRows,faceRows);
reduced.strainOperator=fullModel.strainOperator*Pu;
reduced.bodyPressure=reduced.stiffness\(reduced.gapOperator'*diag(reduced.weight*p.wetAreaFraction));
reduced.bodyTemperature=reduced.stiffness\reduced.thermalLoad;
reduced.bodyTraction=reduced.stiffness\(reduced.gapOperator'*diag(reduced.weight));
reduced.gapCompliance=reduced.gapOperator*reduced.bodyTraction;
reduced.displacementPressure=reduced.bodyPressure; reduced.displacementTemperature=reduced.bodyTemperature;
reduced.symmetry=struct('temperatureLift',Pt,'displacementLift',Pu,'cellLift',jointMap,'multiplicity',8,'fullModel',fullModel);
reduced.scope='One-sector exact cyclic projection for uniform nominal verification; heterogeneous samples require eight sectors.';
pReduced=p;
if ~isscalar(p.aperture0_m)
    field=reshape(p.aperture0_m,nq,8);
    assert(max(abs(field-field(:,1)),[],'all')<1e-15,'PhaseSwitch:Symmetry','Aperture field is not eightfold identical.');
    pReduced.aperture0_m=field(:,1);
end
end
