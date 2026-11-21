function [update,info] = premelting_newton_direction(matrix,rhs,fraction,n,eliminate)
%PREMELTING_NEWTON_DIRECTION Exact elimination of trivial phase equations.
% A projection-inactive phase row contains only its own fraction derivative.
% Solve that equation exactly and remove its row AND column. The remaining
% least-squares objective is identical: removed rows have no free-variable
% derivatives. General mixed-phase rows and additional bound pins retain the
% original rectangular least-squares solve and column scaling.
if nargin<5, eliminate=true; end
nz=numel(rhs); phaseIDs=n+(1:3*n)';
pinned=false(3*n,1); target=zeros(3*n,1); rows=(1:nz)';
if eliminate
    block=matrix(phaseIDs,:);
    diagonal=block(sub2ind(size(block),(1:3*n)',phaseIDs));
    block(sub2ind(size(block),(1:3*n)',phaseIDs))=0;
    trivial=~any(block~=0,2) & diagonal~=0;
    pinned=trivial;
    target(trivial)=rhs(phaseIDs(trivial))./diagonal(trivial);
    rows=rows(~ismember(rows,phaseIDs(trivial)));
end
for attempt=1:3*n+1
    fixed=phaseIDs(pinned); free=setdiff((1:nz)',fixed);
    update=zeros(nz,1); update(fixed)=target(pinned);
    reduced=rhs(rows)-matrix(rows,fixed)*update(fixed);
    % Keep the original column norm, including any removed zero entries.
    columnNorm=sqrt(sum(matrix(:,free).^2,1))';
    assert(all(isfinite(columnNorm)&columnNorm>0),'ThermoInterface:LinearSystem', ...
        'A Newton column has no finite nonzero derivative.');
    update(free)=(matrix(rows,free)./columnNorm'\reduced)./columnNorm;
    phaseUpdate=update(phaseIDs);
    lower=fraction<=1e-13 & phaseUpdate<0;
    upper=fraction>=1-1e-13 & phaseUpdate>0;
    add=(lower|upper)&~pinned;
    if ~any(add), break; end
    pinned=pinned|add; target(add)=double(upper(add))-fraction(add);
end
info=struct('originalSize',nz,'remainingRows',numel(rows), ...
    'remainingColumns',numel(free),'boundIterations',attempt);
end
