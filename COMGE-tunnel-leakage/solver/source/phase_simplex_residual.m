function [r,Jtheta,Jomega] = phase_simplex_residual(theta,omega,scale)
%PHASE_SIMPLEX_RESIDUAL Semismooth equilibrium with three nonnegative areas.
% r=theta-project_simplex(theta-Omega/scale). Subtracting Omega_l is a
% gauge choice only. No phase is mandatory; pure vapour is admissible.
n=size(theta,1); r=zeros(n,3);
if nargout>1, Jtheta=zeros(n,3,3); Jomega=Jtheta; end
assert(isscalar(scale)||numel(scale)==n,'PhaseSwitch:Shape','Supply one positive scale per cell or one scalar.');
assert(all(isfinite(scale)&scale>0),'PhaseSwitch:Domain','Phase scales must be positive and finite.');
if isscalar(scale), scale=repmat(scale,n,1); else, scale=scale(:); end
for j=1:n
    x=theta(j,:)'-(omega(j,:)'-omega(j,1))/scale(j);
    sorted=sort(x,'descend'); cs=cumsum(sorted)-1;
    activeCount=find(sorted-cs./(1:3)'>0,1,'last'); offset=cs(activeCount)/activeCount;
    projected=max(x-offset,0); r(j,:)=theta(j,:)-projected';
    if nargout>1
        active=projected>0; P=zeros(3);
        P(active,active)=eye(sum(active))-ones(sum(active))/sum(active);
        Jtheta(j,:,:)=eye(3)-P; Jomega(j,:,:)=P/scale(j);
    end
end
end
