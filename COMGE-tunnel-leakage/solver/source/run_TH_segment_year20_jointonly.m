function result=run_TH_segment_year20_jointonly(user)
%RUN_TH_SEGMENT_YEAR20_JOINTONLY Current full-history joint-only THM entry.
% Uses v16 and schema-11 annual continuation, including installed U/T.
% Options: model, maxYears, resume, outputFolder. See the THM README.
% For the verified faster numerical preset use run_TH_segment_year20_fast_thm.
% Historical v14 TH checkpoints use run_TH_segment_year20_regime explicitly.
if nargin<1, user=struct(); end
result=run_TH_segment_year20_thm(user);
end
