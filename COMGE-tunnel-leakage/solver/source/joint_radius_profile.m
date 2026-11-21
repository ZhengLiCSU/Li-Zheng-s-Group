function [rq,relativeMobility]=joint_radius_profile(settings,nj,nq)
%JOINT_RADIUS_PROFILE Fixed local geometry with a shared connected pore area.
% Reference radius fixes area/volume with Cref. Local radius changes conductance
% by (r/rref)^2 at fixed connected area, not by independently changing storage.
r=settings.flowWeightedRadius_m;
assert(isnumeric(r)&&isreal(r)&&all(isfinite(r(:)))&&all(r(:)>0) ...
    &&(isscalar(r)||numel(r)==nj),'JointField:Reference','Invalid reference radius.');
r=r(:);if isscalar(r),r=repmat(r,nj,1);end
rq=repmat(r,1,nq);relativeMobility=ones(nj,nq);
if isempty(settings.flowWeightedRadiusProfile_m),return;end
value=settings.flowWeightedRadiusProfile_m;
assert(isnumeric(value)&&isreal(value)&&isequal(size(value),[nj nq]) ...
    &&all(isfinite(value(:)))&&all(value(:)>0),'JointField:Profile', ...
    'Provide one strictly positive radius at each actual joint/path sample.');
rq=double(value);relativeMobility=(rq./r).^2;
assert(all(isfinite(relativeMobility(:)))&&all(relativeMobility(:)>0), ...
    'JointField:Range','Radius mobility overflows or underflows.');
end
