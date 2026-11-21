function report=fit_joint_cycle_coarsening(data,provenance,options)
%FIT_JOINT_CYCLE_COARSENING Fit thawed conductance ratios with diagnostics.
% DATA table: cumulativeDose [1], conductanceRatio=Cref(N)/Cref(0) [1].
% Only zero initial severity, fully thawed observations at fixed compression
% and corrected to one viscosity reference are supported. Dose is NOT years.
% Required provenance: kind, source, material, thawedReferenceCorrected,
% constantCompression, uniformPathDose. Mean dose of heterogeneous path
% points cannot replace their nonlinear severity mean. No automatic
% transfer from rock/concrete to EPDM.
% Profile optimization is toolbox-free; the 2 parameters may be inseparable
% over a short early-cycle record. Fit success is NOT physical validation.
if nargin<3,options=struct();end
o=struct('rateBounds',[1e-5 5],'alphaBounds',[0 5],'gridPoints',121, ...
    'minimumScaledSingularRatio',1e-3,'baselineRatioTolerance',.02);
assert(isstruct(options)&&isscalar(options),'JointCycleFit:Options','Options must be scalar struct.');
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'JointCycleFit:Options','Unknown option: %s.',key{1});
    o.(key{1})=options.(key{1});
end
validateattributes(o.rateBounds,{'numeric'},{'real','finite','positive','increasing','numel',2});
validateattributes(o.alphaBounds,{'numeric'},{'real','finite','nonnegative','increasing','numel',2});
validateattributes(o.gridPoints,{'numeric'},{'scalar','integer','>=',11});
validateattributes(o.minimumScaledSingularRatio,{'numeric'},{'scalar','positive','<',1});
validateattributes(o.baselineRatioTolerance,{'numeric'},{'scalar','positive','<',1});
assert(istable(data)&&all(ismember({'cumulativeDose','conductanceRatio'}, ...
    data.Properties.VariableNames))&&height(data)>=4, ...
    'JointCycleFit:Data','At least four unit-labelled observations are required.');
x=data.cumulativeDose(:);ratio=data.conductanceRatio(:);
validateattributes(x,{'numeric'},{'real','finite','nonnegative'});
validateattributes(ratio,{'numeric'},{'real','finite','positive'});
assert(numel(unique(x(x>0)))>=3&&any(x==0), ...
    'JointCycleFit:Coverage','Require a zero-dose baseline and three distinct positive doses.');
assert(abs(mean(ratio(x==0))-1)<=o.baselineRatioTolerance, ...
    'JointCycleFit:Baseline','Normalize to the same undamaged thawed reference conductance.');
assert(isstruct(provenance)&&isscalar(provenance)&&all(isfield(provenance, ...
    {'kind','source','material','thawedReferenceCorrected','constantCompression','uniformPathDose'})) ...
    &&any(strcmp(provenance.kind,{'experiment','literature_table','literature_digitized','synthetic'})) ...
    &&isscalar(string(provenance.source))&&strlength(string(provenance.source))>0 ...
    &&isscalar(string(provenance.material))&&strlength(string(provenance.material))>0, ...
    'JointCycleFit:Provenance','Explicit data kind, source, material and test conditions are required.');
assert(isequal(provenance.thawedReferenceCorrected,true) ...
    &&isequal(provenance.constantCompression,true)&&isequal(provenance.uniformPathDose,true), ...
    'JointCycleFit:Conditions','Require thawed referenced flow, fixed compression and uniform known path dose.');
z=sqrt(ratio)-1;
grid=linspace(log(o.rateBounds(1)),log(o.rateBounds(2)),o.gridPoints);
profile=zeros(size(grid));alphas=profile;
for k=1:numel(grid),[profile(k),alphas(k)]=objective(grid(k));end
[~,best]=min(profile);lo=grid(max(1,best-1));hi=grid(min(numel(grid),best+1));
[logRate,value]=fminbnd(@objective,lo,hi,optimset('Display','off','TolX',1e-10));
if profile(best)<value,logRate=grid(best);end
[objectiveValue,alpha]=objective(logRate);rate=exp(logRate);
D=joint_cycle_severity(x,rate,0);prediction=(1+alpha*D).^2;
g=-expm1(-rate*x);J=[g,alpha*rate*x.*exp(-rate*x)];
norms=sqrt(sum(J.^2,1));scaled=J./max(norms,realmin);
singular=svd(scaled,0);sratio=singular(end)/max(singular(1),realmin);
boundary=rate<=o.rateBounds(1)*(1+1e-5)||rate>=o.rateBounds(2)*(1-1e-5) ...
    ||alpha<=o.alphaBounds(1)+1e-8||alpha>=o.alphaBounds(2)-1e-8;
identifiable=all(norms>0)&&sratio>=o.minimumScaledSingularRatio&&~boundary;
scopeMatch=strcmp(provenance.material,'gasketed_segment_joint');
report=struct('schemaVersion',1,'provenance',provenance,'options',o, ...
    'doseRate',rate,'radiusGrowthAtUnitDamage',alpha,'initialSeverity',0, ...
    'objectiveInSqrtRatio',objectiveValue,'rmseConductanceRatio',sqrt(mean((prediction-ratio).^2)), ...
    'scaledSingularRatio',sratio,'identifiable',identifiable,'parameterAtBound',boundary, ...
    'materialMatchesSegmentGasket',scopeMatch, ...
    'isSynthetic',strcmp(provenance.kind,'synthetic'),'isIndependentValidation',false, ...
    'profile',table(exp(grid(:)),alphas(:),profile(:), ...
        'VariableNames',{'doseRate','optimalAlpha','objectiveInSqrtRatio'}), ...
    'comparison',table(x,ratio,prediction,prediction-ratio, ...
        'VariableNames',{'cumulativeDose','measuredRatio','predictedRatio','residual'}), ...
    'scope','Fit of the assumed constant-volume law only. Joint geometry/volume assumption needs independent observations; no confidence interval or automatic experimental-validation claim.', ...
    'warnings',{{}});
if ~identifiable,report.warnings{end+1}='Rate and radius mapping are not separately supported by this dose span/bounds.';end
if ~scopeMatch,report.warnings{end+1}='Material differs from a gasketed segment joint; these estimates must not be imported as gasket properties.';end
if report.isSynthetic,report.warnings{end+1}='Synthetic recovery qualifies fitting software only; it provides no literature or experimental calibration.';end

    function [value,a]=objective(logK)
        shape=-expm1(-exp(logK)*x);den=shape'*shape;
        a=min(o.alphaBounds(2),max(o.alphaBounds(1),(shape'*z)/max(den,realmin)));
        value=sum((a*shape-z).^2);
    end
end
