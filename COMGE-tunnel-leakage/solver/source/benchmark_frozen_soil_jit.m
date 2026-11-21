function report=benchmark_frozen_soil_jit(kernel,n,p)
%BENCHMARK_FROZEN_SOIL_JIT Parity and steady-call cost, without solver changes.
% REPORT=benchmark_frozen_soil_jit(KERNEL,N,P) compares the original public
% properties-only law with the isolated numeric JIT entry and a struct adapter.
% The adapter timing includes reconstructing every original property field and
% its empty tangent fields. The numeric-only speedup is NOT the solver speedup.
% Warm-up is separate. Timeit is serial; no parfor or nested OpenMP is used.
% No full TH/THM solve, mesh change, reduced quadrature or code dispatch occurs.

if nargin<1||isempty(kernel), kernel=str2func('frozen_soil_properties_jit_mex'); end
if nargin<2||isempty(n), n=200000; end
if nargin<3||isempty(p), p=fixture(); end
validateattributes(n,{'double'},{'scalar','integer','positive','finite'});
assert(isa(kernel,'function_handle'),'FrozenSoilJIT:Kernel','Supply the isolated entry handle.');
[v,schema]=frozen_soil_jit_parameters(p);
i=(1:n)'; T=-3+9*sin(.017*i); h=-20+28*cos(.013*i);
reference=@() frozen_soil_constitutive(T,h,p,false);
numeric=@() kernel(T,h,v);
adapted=@() asState(kernel(T,h,v),schema);
a=reference(); Y=numeric(); b=adapted();
propertyRelativeError=zeros(1,12);
for j=1:12
    name=schema.propertyNames{j}; expected=a.(name);
    delta=abs(Y(:,j)-expected);
    % Normalize each point with its verification scale. For K the near-zero
    % scale is 2.5e-288, so small frozen K is not hidden by a macroscopic floor.
    denominator=max(abs(expected),propertyAbsoluteTolerance(j)/(4e-13));
    propertyRelativeError(j)=max(delta./denominator);
    assert(all(delta<=4e-13*abs(expected)+propertyAbsoluteTolerance(j)), ...
        'FrozenSoilJIT:Parity','The compiled property %s differs from the original law.',name);
    assert(isequaln(b.(name),Y(:,j)),'FrozenSoilJIT:Adapter','Adapter altered a property.');
end
% Trigger any JIT generation on the first invocation before timing.
reference(); numeric(); adapted();
referenceSeconds=timeit(reference,1);
numericSeconds=timeit(numeric,1);
adaptedSeconds=timeit(adapted,1);
report=struct('matlabRelease',version('-release'),'nPoints',n, ...
    'referenceSeconds',referenceSeconds,'numericJitSeconds',numericSeconds, ...
    'adaptedJitSeconds',adaptedSeconds, ...
    'numericSpeedup',referenceSeconds/numericSeconds, ...
    'adaptedSpeedup',referenceSeconds/adaptedSeconds, ...
    'propertyNames',{schema.propertyNames}, ...
    'propertyVerificationScaledRelativeError',propertyRelativeError, ...
    'kernel',func2str(kernel), ...
    'scope','Pointwise properties only; end-to-end THM speedup is not measured.');
fprintf('%d points: MATLAB %.6g s, JIT numeric %.6g s, JIT adapted %.6g s.\n', ...
    n,referenceSeconds,numericSeconds,adaptedSeconds);
fprintf('Numeric speedup %.3fx; including struct adaptation %.3fx.\n', ...
    report.numericSpeedup,report.adaptedSpeedup);
end

function a=asState(Y,schema)
for j=1:12, a.(schema.propertyNames{j})=Y(:,j); end
derivatives={'slT','slH','siT','siH','dMassT','dMassH', ...
    'dEnergyT','dEnergyH','dkT','dkH','dKT','dKH'};
for j=1:numel(derivatives), a.(derivatives{j})=[]; end
end

function tol=propertyAbsoluteTolerance(j)
if j==12, tol=1e-300; else, tol=2e-14; end
end

function p=fixture()
p=struct('rho',1000,'g',9.81,'L',3.34e5,'Tm',273.15, ...
    'porosity',.35,'thetaR',.05,'vgAlpha',.65,'vgN',1.67, ...
    'Ksat',9.81e-8,'Ss',1e-5,'iceImpedance',3.2, ...
    'cw',4180,'ci',2100,'rhoS',2650,'cs',800, ...
    'kw',.57,'ki',2.2,'ks',2.3,'ka',.025, ...
    'soilViscosity',struct('enabled',true,'referenceTemperature_C',20));
end
