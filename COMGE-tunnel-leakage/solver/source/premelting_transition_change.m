function change = premelting_transition_change(old,next)
%PREMELTING_TRANSITION_CHANGE Monitor physical ice pressure on active domains.
active=old.areaFraction(:,2)>1e-10 | next.areaFraction(:,2)>1e-10;
ice=0;
if any(active), ice=max(abs(next.eos.pressure(active,2)-old.eos.pressure(active,2))); end
change=struct('phase',max(abs(next.areaFraction-old.areaFraction),[],'all'), ...
    'liquidAndWallPressure_Pa',max(abs([next.liquidPressure-old.liquidPressure;next.pressure-old.pressure])), ...
    'activeIcePressure_Pa',ice,'temperature_K',max(abs(next.x(:,4)-old.x(:,4))));
end
