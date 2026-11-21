function proposal=premelting_step_proposal(dt,phaseChange,pressureChange,damageChange,icePressureChange,o)
%PREMELTING_STEP_PROPOSAL Predict next trial; never bypass acceptance gates.
% Rates are inferred only from the committed step. Abrupt transitions still
% undergo the original reject/halve loop. Maximum step and all physical,
% conservation and transition limits are unchanged.
factor=1.5;
if strcmp(o.stepController,'transition_predictor')
    use=max([phaseChange/o.maximumPhaseIncrement, ...
        pressureChange/o.maximumPressureIncrement_Pa, ...
        damageChange/o.maximumDamageIncrement, ...
        icePressureChange/o.maximumActiveIcePressureIncrement_Pa]);
    if use>0, factor=min(factor,.9/use); end
end
trial=factor*dt;
if isfield(o,'minimumStepSeconds'), trial=max(o.minimumStepSeconds,trial); end
proposal=min(o.dtDays*86400,trial);
end
