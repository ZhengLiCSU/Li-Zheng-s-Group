function a=scale_soil_conductivity_state(a,multiplier)
%SCALE_SOIL_CONDUCTIVITY_STATE Scale K and ALL its available exact partials.
% A fixed material multiplier changes hydraulic mobility, not phase storage,
% retention or thermal conductivity. Application follows the THM bulk law.
names={'K','dKT','dKH','dKEps'};
for j=1:numel(names)
    if isfield(a,names{j})&&~isempty(a.(names{j}))
        a.(names{j})=a.(names{j}).*multiplier;
    end
end
end
