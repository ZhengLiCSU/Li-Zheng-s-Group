function [outer,inner,head] = premelting_boundary_at(time,referenceT,o)
%PREMELTING_BOUNDARY_AT Synthetic cycle or a supplied SI time series.
if isempty(o.boundarySchedule)
    outer=referenceT;
    inner=referenceT-o.temperatureAmplitude_K*sin(pi*time/(o.periodDays*86400))^2;
    head=o.reservoirPressure_Pa;
else
    s=o.boundarySchedule;
    assert(time>=s.time_s(1) && time<=s.time_s(end),'PhaseSwitch:Input','Boundary extrapolation is forbidden.');
    outer=interp1(s.time_s,s.outerTemperature_K,time,'linear');
    inner=interp1(s.time_s,s.innerTemperature_K,time,'linear');
    head=interp1(s.time_s,s.reservoirPressure_Pa,time,'linear');
end
end
