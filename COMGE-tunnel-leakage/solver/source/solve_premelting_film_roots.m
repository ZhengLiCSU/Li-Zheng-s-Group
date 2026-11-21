function [root,info] = solve_premelting_film_roots(split,low,high,f)
%SOLVE_PREMELTING_FILM_ROOTS Batched safeguarded Newton on the stable branch.
% Same Pi(h)=split, positive bounded domain and 1e-22 m absolute root
% tolerance as the scalar reference. No interpolated/approximated potential.
% Always retain a sign bracket; fall back to bisection for unsafe Newton
% steps. This helper only solves interior candidates, not their global
% energy selection or accepted-state thin-film-domain audit.
arguments
    split (:,1) double {mustBeFinite}
    low (:,1) double {mustBePositive,mustBeFinite}
    high (:,1) double {mustBePositive,mustBeFinite}
    f (1,1) struct
end
if isscalar(low), low=repmat(low,size(split)); end
if isscalar(high), high=repmat(high,size(split)); end
assert(isequal(size(low),size(split),size(high)) && all(high>low), ...
    'Premelting:RootBracket','Root bounds must match and increase.');
[fl,~]=balance(low,split,f); [fh,~]=balance(high,split,f);
assert(all(fl>0 & fh<0),'Premelting:RootBracket', ...
    'Stable interior roots require strictly opposite endpoint signs.');
root=(low+high)/2;
positive=split>0 & f.shortRangeEnergy_Jm2>0;
guess=-log(split(positive)/(f.inverseLength_minv*f.shortRangeEnergy_Jm2))/f.inverseLength_minv;
root(positive)=min(high(positive),max(low(positive),guess));
active=true(size(root)); tolerance=1e-22; count=zeros(size(root));
for iteration=1:160
    ids=find(active); if isempty(ids), break; end
    x=root(ids); [value,slope]=balance(x,split(ids),f); count(ids)=count(ids)+1;
    zero=value==0; active(ids(zero))=false;
    lower=value>0; upper=value<0;
    low(ids(lower))=x(lower); fl(ids(lower))=value(lower);
    high(ids(upper))=x(upper); fh(ids(upper))=value(upper);
    narrow=high(ids)-low(ids)<=tolerance;
    done=~zero & narrow;
    best=low(ids(done)); useUpper=abs(fh(ids(done)))<abs(fl(ids(done)));
    selected=ids(done); best(useUpper)=high(selected(useUpper));
    root(selected)=best; active(selected)=false;
    keep=~zero & ~narrow; ids=ids(keep); x=x(keep); value=value(keep); slope=slope(keep);
    if isempty(ids), continue; end
    width=high(ids)-low(ids); proposal=x-value./slope;
    safe=isfinite(proposal) & slope<0 & ...
        proposal>low(ids)+.05*width & proposal<high(ids)-.05*width;
    midpoint=(low(ids)+high(ids))/2; proposal(~safe)=midpoint(~safe);
    root(ids)=proposal;
end
assert(~any(active),'Premelting:RootConvergence','Batched stable-root bracket did not converge.');
[defect,~]=balance(root,split,f);
info=struct('iterations',count,'finalBracketWidth_m',high-low,'pressureDefect_Pa',defect);
end

function [value,slope]=balance(h,split,f)
e=f.shortRangeEnergy_Jm2*exp(-f.inverseLength_minv*h);
value=f.inverseLength_minv*e-f.hamaker_J./(6*pi*h.^3)-split;
slope=-f.inverseLength_minv^2*e+f.hamaker_J./(2*pi*h.^4);
end
