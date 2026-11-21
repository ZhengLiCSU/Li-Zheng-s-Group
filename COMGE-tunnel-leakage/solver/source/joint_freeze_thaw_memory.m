function [s,event]=joint_freeze_thaw_memory(s,ice,timeDay,active,c)
%JOINT_FREEZE_THAW_MEMORY Accepted-state local liquid/ice cycle tracker.
% Each path sample must first thaw, then freeze, then thaw again. A cold
% initial half-cycle is ignored. No calendar cycles or rejected iterations.
% Completion dose = ((peakIce-thawFraction)/(1-thawFraction))^exponent.
% Dsample=1-(1-Dinitial)*exp(-doseRate*cumulativeDose); joint D is path mean.
% Partial cycles, duration and peaks persist across annual checkpoints.
c=joint_cycle_damage_defaults(c);
assert(isnumeric(ice)&&isreal(ice)&&ismatrix(ice)&&all(isfinite(ice(:))) ...
    &&all(ice(:)>=0&ice(:)<=1)&&isscalar(timeDay)&&isfinite(timeDay), ...
    'JointCycle:State','Supply finite sample ice fractions in [0,1] and time.');
assert((isnumeric(active)||islogical(active))&&isreal(active) ...
    &&all(isfinite(active(:)))&&all(ismember(active(:),[0 1])), ...
    'JointCycle:Mask','Active mask must contain only finite 0/1 values.');
[nj,nq]=size(ice);active=logical(active(:));
if isscalar(active)&&nj>1,active=repmat(active,nj,1);end
assert(numel(active)==nj&&nq>0,'JointCycle:Shape','Active mask/path dimensions differ.');
rate=expand(c.doseRate,nj);D0=expand(c.initialSeverity,nj);
event=struct('completed',false(nj,nq),'dose',zeros(nj,nq));
if isempty(s)
    s=struct('revision','local_ice_cycle_v1_20261002','lastDay',timeDay, ...
        'armed',ice<=c.thawFraction&active,'frozen',false(nj,nq), ...
        'peakIce',zeros(nj,nq),'freezeStartDay',nan(nj,nq), ...
        'cycleCount',zeros(nj,nq),'cumulativeDose',zeros(nj,nq), ...
        'sampleSeverity',repmat(D0,1,nq),'targetSeverity',D0, ...
        'lastIce',ice);
    return
end
assert(isstruct(s)&&strcmp(s.revision,'local_ice_cycle_v1_20261002') ...
    &&isequal(size(s.lastIce),size(ice))&&timeDay>s.lastDay, ...
    'JointCycle:Memory','Memory dimensions/revision/time are inconsistent.');
validate_joint_cycle_memory(s,nj,nq,c);
% Completed cycles use the previously armed freezing excursion. Warm and
% cold chatter cannot generate repeated cycles within one accepted sample.
enter=s.armed&~s.frozen&ice>=c.freezeFraction&active;
s.frozen(enter)=true;s.armed(enter)=false;s.freezeStartDay(enter)=timeDay;
s.peakIce(s.frozen)=max(s.peakIce(s.frozen),ice(s.frozen));
close=s.frozen&ice<=c.thawFraction&active;
valid=close&(timeDay-s.freezeStartDay>=c.minimumFrozenDays);
event.completed=valid;
event.dose(valid)=((s.peakIce(valid)-c.thawFraction)/(1-c.thawFraction)).^c.doseExponent;
s.cycleCount(valid)=s.cycleCount(valid)+1;
s.cumulativeDose=s.cumulativeDose+event.dose;
s.frozen(close)=false;s.peakIce(close)=0;s.freezeStartDay(close)=NaN;
s.armed=s.armed|(ice<=c.thawFraction&~s.frozen&active);
s.sampleSeverity=joint_cycle_severity(s.cumulativeDose,rate,D0);
s.targetSeverity=mean(s.sampleSeverity,2);
s.lastIce=ice;s.lastDay=timeDay;
end
function v=expand(v,n)
assert(isscalar(v)||numel(v)==n,'JointCycle:Shape','Settings require scalar or joint vector.');
v=double(v(:));if isscalar(v),v=repmat(v,n,1);end
end
