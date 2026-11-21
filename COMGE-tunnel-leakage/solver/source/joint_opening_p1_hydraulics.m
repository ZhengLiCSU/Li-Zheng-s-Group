function [C,DCU,DCT,D]=joint_opening_p1_hydraulics(Cref,B,ur,Tprofile,options)
%JOINT_OPENING_P1_HYDRAULICS Continuous endpoint-aware P1 slit resistance.
% [C,DCU,DCT,D]=joint_opening_p1_hydraulics(Cref,Bgap,ur,Tprofile,options)
% uses the ACTUAL radial mechanical nodes in build_joint_gap_profiles.
% The new slit aperture is b(r)=f*max(g_P1(r),0), only where Cref>0.
% An open path integrates cubic resistance EXACTLY in aperture, splitting
% mechanical intervals at thermal-cell boundaries. Temperature/mobility is
% piecewise constant in each thermal length cell; it is not interpolated
% into the unresolved endpoint. For positive endpoint apertures a,b,
%       integral dr/b(r)^3 = length*(a+b)/(2*a^2*b^2).
% Thus Copen tends continuously to zero when any endpoint closes, unlike a
% midpoint-only aperture rule followed by an independent nodal-closure mask.
% Genuine closure/zero thermal mobility blocks only the NEW flow branch.
% Partial-opening geometric inventory remains positive and is integrated
% exactly. No aperture/conductance floor or geometric smoothing is used.
%
% C is nj-by-1; DCU is sparse nj-by-numel(ur); DCT is nj-by-nq, the partial
% with respect to thermal SAMPLES. With Pthermal joint-major interpolation,
% global CT = sparse(repelem((1:nj)',nq),(1:nj*nq)', ...
%             reshape(DCT',[],1),nj,nj*nq)*Pthermal.
% D.dC_dNodalGap is nj-by-nr in the ORIGINAL mechanical radial-node order.
% Cref*baseMobilityFactor is a separate capillary branch; its optional
% baseMobilityDerivative_dT is nj-by-nq, exactly as in the old profile API.
%
% The minimum exact-resistance version requires wallEnergyDifference_J_m2=0
% (the scalar helper's default). A nonzero aperture-dependent Gibbs-Thomson
% temperature requires spatial thermodynamic integration and is REJECTED,
% not silently averaged or approximated by this kernel. Other constitutive
% settings are validated/reused by joint_opening_hydraulics. The old scalar
% and midpoint-profile helpers remain unchanged and backwards compatible.
%
% Q=joint_opening_p1_hydraulics('build',Bgap) constructs static sparse piece
% maps. build_joint_gap_profiles caches Q as B.resistanceIntegration; callers
% with an older cache are supported by rebuilding these small maps once per
% call. Thermal cells have lengths L*B.sampleWeights, outer-to-inner.
% Returned volumes are diagnostics for the existing zero-storage limit,
% not a moving-channel continuity or pressure/storage equation.

if nargin==2 && (ischar(Cref)||(isstring(Cref)&&isscalar(Cref)))
    assert(strcmp(Cref,'build'),'JointOpeningP1:Command','Unknown command.');
    C=buildPieces(B); DCU=[]; DCT=[]; D=[]; return
end
if nargin<4, Tprofile=[]; end
if nargin<5 || isempty(options), options=struct(); end
assert(isstruct(options)&&isscalar(options),'JointOpeningP1:Options', ...
    'Options must be a scalar struct.');
if isfield(B,'resistanceIntegration')
    Q=B.resistanceIntegration;
else
    Q=buildPieces(B);
end
nj=Q.nj; nr=Q.nr; nq=Q.nq; np=nj*nq;
assert(isnumeric(ur)&&isreal(ur)&&isvector(ur) ...
    &&numel(ur)==size(B.nodalGapU,2)&&all(isfinite(ur(:))), ...
    'JointOpeningP1:Displacement','UR must contain all finite ring DOFs.');
ur=double(ur(:));
Cref=jointValue(Cref,nj,0,Inf,'Conductance');
base=ones(nj,1); fraction=zeros(nj,1);
if isfield(options,'baseMobilityFactor')
    base=jointValue(options.baseMobilityFactor,nj,0,Inf,'BaseFactor');
end
if isfield(options,'openingTransmissionFraction')
    fraction=jointValue(options.openingTransmissionFraction,nj,0,1,'Fraction');
end
baseT=zeros(nj,nq);
if isfield(options,'baseMobilityDerivative_dT') && ~isempty(options.baseMobilityDerivative_dT)
    baseT=options.baseMobilityDerivative_dT;
    assert(isnumeric(baseT)&&isreal(baseT)&&isequal(size(baseT),[nj nq]) ...
        &&all(isfinite(baseT(:))),'JointOpeningP1:BaseDerivative', ...
        'The base-factor derivative must be nj-by-nq.');
    baseT=double(baseT);
end
mask=true(nj,1);
if isfield(options,'openPathMask')
    mask=jointValue(double(options.openPathMask),nj,0,1,'OpenPathMask');
    assert(all(ismember(mask,[0 1])),'JointOpeningP1:OpenPathMask','Mask must be logical 0/1.');
    mask=logical(mask);
end
thermalProvided=~isempty(Tprofile);
if ~thermalProvided
    frozen=isfield(options,'freezeEnabled')&&options.freezeEnabled;
    viscous=isfield(options,'viscosityEnabled')&&options.viscosityEnabled;
    assert(~frozen&&~viscous,'JointOpeningP1:Temperature', ...
        'Enabled freezing or viscosity requires temperatures.');
    Tprofile=zeros(nj,nq);
end
assert(isnumeric(Tprofile)&&isreal(Tprofile)&&isequal(size(Tprofile),[nj nq]) ...
    &&all(isfinite(Tprofile(:))),'JointOpeningP1:Temperature', ...
    'Temperatures must be finite nj-by-nq samples.');
if isfield(options,'wallEnergyDifference_J_m2')
    gamma=options.wallEnergyDifference_J_m2;
    assert(isnumeric(gamma)&&isreal(gamma)&&isscalar(gamma)&&isfinite(gamma)&&gamma==0, ...
        'JointOpeningP1:ApertureDependentFreezing', ...
        'Exact P1 resistance currently requires zero wall-energy difference.');
end
if isfield(options,'pathLength_m')
    supplied=jointValue(options.pathLength_m,nj,realmin,Inf,'PathLength');
    assert(all(abs(supplied-Q.pathLength_m)<=1e-10*max(1,Q.pathLength_m)), ...
        'JointOpeningP1:PathLength','Constitutive and actual mechanical lengths disagree.');
end
if isfield(options,'sampleWeights') && ~isempty(options.sampleWeights)
    supplied=double(options.sampleWeights(:)');
    assert(numel(supplied)==nq&&all(isfinite(supplied))&&all(supplied>0), ...
        'JointOpeningP1:Weights','Invalid thermal length fractions.');
    supplied=supplied/sum(supplied);
    assert(max(abs(supplied-Q.sampleWeights))<1e-12,'JointOpeningP1:Weights', ...
        'Thermal length fractions must match the cached geometry.');
end
% Mobility is independent of aperture when wall-energy difference is zero.
% Reuse the validated scalar constitutive law on a positive unit aperture;
% dummy conductance only enables evaluation and never enters the result.
localOptions=options;
if isfield(localOptions,'openPathMask'), localOptions=rmfield(localOptions,'openPathMask'); end
localOptions.openingTransmissionFraction=1;
localOptions.pathLength_m=repelem(Q.pathLength_m,nq);
localOptions.sampleWeights=1; localOptions.baseMobilityFactor=0;
localOptions.baseMobilityDerivative_dT=zeros(np,1);
[~,~,localDT,localD]=joint_opening_hydraulics(ones(np,1),ones(np,1), ...
    reshape(double(Tprofile)',[],1),localOptions);
mobility=localD.localMobility(:);
mobilityT=localDT./localD.referenceOpenConductance;
phase=localD.liquidFractionBySample(:);
phaseT=localD.dLiquidFraction_dTemperature(:);
s=localD.settings;
muReference=water_viscosity_iapws(s.referenceTemperature_C);
logCoefficient=log(s.waterDensity_kg_m3)+log(s.gravity_m_s2)-log(12)-log(muReference);
nodal=double(B.nodalGap0(:))+B.nodalGapU*ur;
assert(all(isfinite(nodal)),'JointOpeningP1:Displacement','Nodal openings overflow.');
nodalJoint=reshape(nodal,nr,nj)';
wholeOpen=all(nodalJoint>0,2);
ga=Q.gapLeftMap*nodal; gb=Q.gapRightMap*nodal;
pieceFraction=fraction(Q.jointID).*double(Cref(Q.jointID)>0);
[integral,integralA,integralB]=positiveP1Integral(ga,gb,Q.length_m);
pieceVolume=pieceFraction.*integral;
volumeGap=spdiags(pieceFraction.*integralA,0,Q.npiece,Q.npiece)*Q.gapLeftMap ...
    +spdiags(pieceFraction.*integralB,0,Q.npiece,Q.npiece)*Q.gapRightMap;
volumeJac=Q.jointPieceMap*volumeGap;
V=Q.jointPieceMap*pieceVolume;
sampleVolume=Q.samplePieceMap*pieceVolume;
openC=zeros(nj,1); referenceOpen=zeros(nj,1);
leftSlope=zeros(Q.npiece,1); rightSlope=leftSlope;
temperatureSlope=leftSlope; logResistance=inf(Q.npiece,1);
pieceResistanceWeight=zeros(Q.npiece,1);
for joint=1:nj
    indices=find(Q.jointID==joint);
    if Cref(joint)==0 || fraction(joint)==0 || ~wholeOpen(joint), continue; end
    thermalIDs=Q.thermalSampleID(indices);
    a=fraction(joint)*ga(indices); b=fraction(joint)*gb(indices);
    la=log(a); lb=log(b); largest=max(la,lb);
    logSum=largest+log1p(exp(min(la,lb)-largest));
    geometryLog=log(Q.length_m(indices))+logSum-log(2)-2*la-2*lb;
    largest=max(geometryLog);
    referenceOpen(joint)=exp(logCoefficient-largest-log(sum(exp(geometryLog-largest))));
    if ~mask(joint) || any(mobility(thermalIDs)<=0), continue; end
    lr=geometryLog-log(mobility(thermalIDs));
    largest=max(lr); normalized=exp(lr-largest);
    logTotal=largest+log(sum(normalized));
    logC=logCoefficient-logTotal;
    openC(joint)=exp(logC);
    logWeight=lr-logTotal;
    scale=max(ga(indices),gb(indices));
    normalizedA=ga(indices)./scale; normalizedB=gb(indices)./scale;
    ratioA=normalizedA./(normalizedA+normalizedB);
    ratioB=normalizedB./(normalizedA+normalizedB);
    leftSlope(indices)=exp(logC+logWeight-log(ga(indices))+log(2-ratioA));
    rightSlope(indices)=exp(logC+logWeight-log(gb(indices))+log(2-ratioB));
    mt=mobilityT(thermalIDs); nonzero=mt~=0;
    tp=indices(nonzero); ids=thermalIDs(nonzero);
    temperatureSlope(tp)=sign(mt(nonzero)).*exp(logC+logWeight(nonzero) ...
        +log(abs(mt(nonzero)))-log(mobility(ids)));
    logResistance(indices)=lr;
    pieceResistanceWeight(indices)=exp(logWeight);
end
nodalDerivative=Q.jointPieceMap*( ...
    spdiags(leftSlope,0,Q.npiece,Q.npiece)*Q.gapLeftMap ...
    +spdiags(rightSlope,0,Q.npiece,Q.npiece)*Q.gapRightMap);
DCU=nodalDerivative*sparse(B.nodalGapU);
sampleDerivative=sparse(Q.jointID,Q.thermalSampleID,temperatureSlope,nj,np);
DCT=reshape(full(sum(sampleDerivative,1)),nq,nj)'+Cref.*baseT;
C=Cref.*base+openC;
openFactor=zeros(nj,1); positiveReference=referenceOpen>0;
openFactor(positiveReference)=openC(positiveReference)./referenceOpen(positiveReference);
status=repmat({'closed_or_disabled'},nj,1);
status(Cref==0)={'sealed'}; status(Cref>0 & wholeOpen & mask)={'ice_blocked_or_underflow'};
status(openC>0)={'open_p1_slit'};
nodalGapDerivative=zeros(nj,nr);
for joint=1:nj
    nodalGapDerivative(joint,:)=full(nodalDerivative(joint,(joint-1)*nr+(1:nr)));
end
liquidPiece=pieceVolume.*phase(Q.thermalSampleID);
icePiece=pieceVolume-liquidPiece;
sampleLiquid=Q.samplePieceMap*liquidPiece;
sampleIce=Q.samplePieceMap*icePiece;
sampleLiquidT=sampleVolume.*phaseT;
D=struct('conductance',C,'baseConductance',Cref.*base,'openConductance',openC, ...
    'referenceOpenConductance',referenceOpen,'openMobilityFactor',openFactor, ...
    'dC_dNodalGap',nodalGapDerivative,'nodalGapTangent',nodalDerivative, ...
    'dC_dU',DCU,'dC_dTemperatureProfile',DCT, ...
    'nodalGap_m',nodalJoint,'nodalAperture_m',fraction.*max(nodalJoint,0).*double(Cref>0), ...
    'wholePathOpen',wholeOpen,'openPathMask',mask,'sealedMask',Cref==0,'status',{status}, ...
    'newslitGeometricVolume_m2',V,'newslitWaterVolume_m2',V, ...
    'dVolume_dU',volumeJac*sparse(B.nodalGapU),'nodalVolumeTangent',volumeJac, ...
    'sampleVolume_m2',reshape(sampleVolume,nq,nj)', ...
    'liquidVolumeBySample_m2',reshape(sampleLiquid,nq,nj)', ...
    'iceVolumeBySample_m2',reshape(sampleIce,nq,nj)', ...
    'dLiquidVolume_dTemperature_m2_K',reshape(sampleLiquidT,nq,nj)', ...
    'dIceVolume_dTemperature_m2_K',-reshape(sampleLiquidT,nq,nj)', ...
    'liquidVolume_m2',Q.jointPieceMap*liquidPiece,'iceVolume_m2',Q.jointPieceMap*icePiece, ...
    'dLiquidVolume_dU',Q.jointPieceMap*spdiags(phase(Q.thermalSampleID),0,Q.npiece,Q.npiece) ...
        *volumeGap*sparse(B.nodalGapU), ...
    'liquidFractionBySample',reshape(phase,nq,nj)', ...
    'localMobility',reshape(mobility,nq,nj)', ...
    'freezeTemperature_C',reshape(localD.freezeTemperature_C,nq,nj)', ...
    'pieceLogResistance',logResistance,'pieceResistanceWeight',pieceResistanceWeight, ...
    'integration',Q,'thermalProfileProvided',thermalProvided,'settings',s, ...
    'apertureIntegration','Exact cubic resistance of actual positive P1 radial apertures; no midpoint aperture surrogate', ...
    'thermalIntegration','Piecewise constant mobility on actual thermal length cells; zero wall-energy difference required', ...
    'storageInterpretation','Exact geometric diagnostic, not a dynamic pressure/inventory equation');
end

function Q=buildPieces(B)
assert(isstruct(B)&&isscalar(B)&&all(isfield(B, ...
    {'nr','nj','nq','nodalRadius_m','nodalGap0','nodalGapU','radiusOuter_m','radiusInner_m','sampleWeights'})), ...
    'JointOpeningP1:Cache','Supply a build_joint_gap_profiles cache.');
nr=B.nr; nj=B.nj; nq=B.nq;
assert(all(isfinite([nr nj nq]))&&nr>=2&&nj>=1&&nq>=1 ...
    &&all([nr nj nq]==fix([nr nj nq])),'JointOpeningP1:Cache','Invalid inventories.');
radius=double(B.nodalRadius_m); outer=double(B.radiusOuter_m(:)); inner=double(B.radiusInner_m(:));
assert(isequal(size(radius),[nr nj])&&all(isfinite(radius(:)))&&all(radius(:)>0) ...
    &&numel(outer)==nj&&numel(inner)==nj&&all(isfinite(outer))&&all(isfinite(inner)) ...
    &&all(outer>inner&inner>0),'JointOpeningP1:Cache','Invalid actual radii.');
assert(numel(B.nodalGap0)==nj*nr&&size(B.nodalGapU,1)==nj*nr, ...
    'JointOpeningP1:Cache','Nodal gap maps must use joint-major rows.');
weights=double(B.sampleWeights(:)');
assert(numel(weights)==nq&&all(isfinite(weights))&&all(weights>0), ...
    'JointOpeningP1:Weights','Thermal lengths must have positive fractions.');
weights=weights/sum(weights); length=outer-inner;
boundary=outer'-[0;cumsum(weights(:))]*length'; boundary(end,:)=inner';
assert(all(diff(boundary,1,1)<0,'all'),'JointOpeningP1:Weights', ...
    'Thermal fractions collapse a cell at floating-point precision.');
maximum=nj*(nr-1)*nq;
leftRow=zeros(2*maximum,1); leftCol=leftRow; leftVal=leftRow;
rightRow=leftRow; rightCol=leftRow; rightVal=leftRow;
jointID=zeros(maximum,1); thermalID=jointID; pieceLength=jointID;
startRadius=jointID; endRadius=jointID; cursor=0;
for joint=1:nj
    [r,order]=sort(radius(:,joint));
    assert(all(diff(r)>0)&&abs(r(1)-inner(joint))<=1e-9*max(1,outer(joint)) ...
        &&abs(r(end)-outer(joint))<=1e-9*max(1,outer(joint)), ...
        'JointOpeningP1:Cache','Mechanical radii and thermal path boundaries disagree.');
    for k=1:nr-1
        for cellID=1:nq
            low=max(r(k),boundary(cellID+1,joint));
            high=min(r(k+1),boundary(cellID,joint));
            if high<=low, continue; end
            cursor=cursor+1; slots=2*(cursor-1)+(1:2);
            columns=(joint-1)*nr+order(k+(0:1));
            al=(low-r(k))/(r(k+1)-r(k)); ar=(high-r(k))/(r(k+1)-r(k));
            leftRow(slots)=cursor; leftCol(slots)=columns; leftVal(slots)=[1-al;al];
            rightRow(slots)=cursor; rightCol(slots)=columns; rightVal(slots)=[1-ar;ar];
            jointID(cursor)=joint; thermalID(cursor)=(joint-1)*nq+cellID;
            pieceLength(cursor)=high-low; startRadius(cursor)=low; endRadius(cursor)=high;
        end
    end
end
slots=1:2*cursor; pieces=1:cursor;
A=sparse(leftRow(slots),leftCol(slots),leftVal(slots),cursor,nj*nr);
Z=sparse(rightRow(slots),rightCol(slots),rightVal(slots),cursor,nj*nr);
jointID=jointID(pieces); thermalID=thermalID(pieces); pieceLength=pieceLength(pieces);
jointMap=sparse(jointID,pieces,1,nj,cursor);
sampleMap=sparse(thermalID,pieces,1,nj*nq,cursor);
assert(max(abs(jointMap*pieceLength-length))<1e-10*max(1,max(length)), ...
    'JointOpeningP1:Cache','Radial resistance pieces must cover every complete joint.');
Q=struct('gapLeftMap',A,'gapRightMap',Z,'length_m',pieceLength, ...
    'jointID',jointID,'thermalSampleID',thermalID,'jointPieceMap',jointMap, ...
    'samplePieceMap',sampleMap,'startRadius_m',startRadius(pieces),'endRadius_m',endRadius(pieces), ...
    'thermalBoundaryRadius_m',boundary,'thermalCellLength_m',-diff(boundary,1,1), ...
    'pathLength_m',length,'sampleWeights',weights,'nj',nj,'nr',nr,'nq',nq,'npiece',cursor);
end

function [I,IA,IB]=positiveP1Integral(a,b,length)
I=zeros(size(a)); IA=I; IB=I;
positive=a>=0 & b>=0 & (a>0|b>0);
I(positive)=length(positive).*(.5*a(positive)+.5*b(positive));
IA(positive)=.5*length(positive); IB(positive)=IA(positive);
left=a>0 & b<0; right=b>0 & a<0;
scale=max(a(left),-b(left)); ratio=(a(left)./scale)./(a(left)./scale-b(left)./scale);
I(left)=.5*length(left).*a(left).*ratio;
IA(left)=.5*length(left).*ratio.*(2-ratio); IB(left)=.5*length(left).*ratio.^2;
scale=max(b(right),-a(right)); ratio=(b(right)./scale)./(b(right)./scale-a(right)./scale);
I(right)=.5*length(right).*b(right).*ratio;
IB(right)=.5*length(right).*ratio.*(2-ratio); IA(right)=.5*length(right).*ratio.^2;
end

function value=jointValue(value,nj,minimum,maximum,name)
assert(isnumeric(value)&&isreal(value)&&isvector(value)&&all(isfinite(value(:))) ...
    &&all(value(:)>=minimum&value(:)<=maximum)&&(isscalar(value)||numel(value)==nj), ...
    ['JointOpeningP1:' name],'Invalid scalar or per-joint %s.',name);
value=double(value(:)); if isscalar(value), value=repmat(value,nj,1); end
end
