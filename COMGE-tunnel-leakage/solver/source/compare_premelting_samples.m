function comparison = compare_premelting_samples(coarse,fine)
%COMPARE_PREMELTING_SAMPLES Fixed 1% target, never replaces independent audit.
keys={'leakage_Lpm','timeIntegratedIceVolume_m2day','timeIntegratedVaporVolume_m2day','peakActiveIcePressure_Pa'};
comparison=struct('relativeToFine',struct(),'tolerance',0.01,'passed',true);
for k=1:numel(keys)
    a=coarse.(keys{k}); b=fine.(keys{k});
    if isnan(a) && isnan(b)
        errorValue=0; % No ice in either trajectory.
    elseif ~isfinite(a) || ~isfinite(b)
        errorValue=Inf;
    else
        errorValue=abs(a-b)/max(abs(b),realmin);
    end
    comparison.relativeToFine.(keys{k})=errorValue;
    comparison.passed=comparison.passed && errorValue<=comparison.tolerance;
end
comparison.scope='Two trajectories for stated meshes/time controls; does not prove convergence across an untested random-field parameter range.';
end
