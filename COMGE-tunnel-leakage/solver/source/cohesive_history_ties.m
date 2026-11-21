function tied = cohesive_history_ties(gap,history,c)
%COHESIVE_HISTORY_TIES Roundoff neighbourhood of a loading/unloading kink.
% Scale is a length, and the factor covers the condensed mechanics solve.
tol=128*eps*max(max(abs(gap),abs(history)),c.onsetGap_m);
tied=c.enabled & history>c.onsetGap_m & history<c.failureGap_m & abs(gap-history)<=tol;
end
