function leak=summarize_segment_leakage(history,leak,startDay)
%SUMMARIZE_SEGMENT_LEAKAGE Integrate accepted-step leakage without resampling bias.
% All volumes are m^3 per metre of tunnel axis; rates are m^2/s.
assert(size(history,1)==size(leak.segmentRate_m2_per_s,1));
assert(size(history,1)==size(leak.jointRate_m2_per_s,1));
assert(size(history,1)==size(leak.bodyRate_m2_per_s,1));
assert(~isempty(history),'At least one accepted step is required.');
assert(all(diff(history(:,1))>0)&&all(history(:,2)>0));
t=[startDay;history(:,1)]; dt=history(:,2)*86400;
assert(max(abs(diff(t)-history(:,2)))<1e-6,'Leakage history has a time gap.');
seg=leak.segmentRate_m2_per_s; joint=leak.jointRate_m2_per_s;
body=leak.bodyRate_m2_per_s;
assert(all(seg(:)>=-1e-14)&&all(joint(:)>=-1e-14)&&all(body(:)>=-1e-14));
assert(max(abs(sum(seg,2)-history(:,10)))<1e-12, ...
    'Segment rates do not sum to the conservative total leakage rate.');
assert(max(abs(body+sum(joint,2)-history(:,10)))<1e-12, ...
    'Body and joint rates do not sum to the conservative total leakage rate.');
leak.totalRate_m2_per_s=sum(seg,2);
leak.volumeBySegment_m3_per_m=sum(seg.*dt,1);
leak.volumeByJoint_m3_per_m=sum(joint.*dt,1);
leak.bodyVolume_m3_per_m=sum(body.*dt);
leak.totalVolume_m3_per_m=sum(history(:,10).*dt);
assert(abs(sum(leak.volumeBySegment_m3_per_m)-leak.totalVolume_m3_per_m)<1e-8);
span=t(end)-startDay;
edges=startDay+(0:floor(span))';
if edges(end)<t(end)-1e-8, edges=[edges;t(end)]; end
if abs(edges(end)-t(end))<1e-8, edges(end)=t(end); end
leak.dailyBinEndDay=edges(2:end);
leak.dailyVolumeBySegment_m3_per_m=integratedBins(t,dt,seg,edges);
leak.dailyVolumeByJoint_m3_per_m=integratedBins(t,dt,joint,edges);
leak.dailyBodyVolume_m3_per_m=integratedBins(t,dt,body,edges);
leak.dailyTotalVolume_m3_per_m=sum(leak.dailyVolumeBySegment_m3_per_m,2);
% Calendar months in a repeated 365-day year, with February fixed at 28 days.
if abs(span-365)<1e-7 && abs(startDay/365-round(startDay/365))<1e-8
    monthDays=[31 28 31 30 31 30 31 31 30 31 30 31];
    monthEdges=[0 cumsum(monthDays)];
    leak.monthlyVolumeBySegment_m3_per_m=zeros(12,size(seg,2));
    leak.monthlyVolumeByJoint_m3_per_m=zeros(12,size(joint,2));
    for month=1:12
        ix=(monthEdges(month)+1):monthEdges(month+1);
        leak.monthlyVolumeBySegment_m3_per_m(month,:)=sum(leak.dailyVolumeBySegment_m3_per_m(ix,:),1);
        leak.monthlyVolumeByJoint_m3_per_m(month,:)=sum(leak.dailyVolumeByJoint_m3_per_m(ix,:),1);
    end
    leak.monthlyTotalVolume_m3_per_m=sum(leak.monthlyVolumeBySegment_m3_per_m,2);
end
end

function bins=integratedBins(t,dt,rate,edges)
if isempty(rate), bins=zeros(numel(edges)-1,0); return; end
cumulative=[zeros(1,size(rate,2));cumsum(rate.*dt,1)];
bins=diff(interp1(t,cumulative,edges,'linear'),1,1);
end
