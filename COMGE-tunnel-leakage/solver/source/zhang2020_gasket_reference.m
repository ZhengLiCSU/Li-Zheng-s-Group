function result = zhang2020_gasket_reference(ageYears, openingMillimetres)
%ZHANG2020_GASKET_REFERENCE Reproduce published 20 C EPDM aging relations.
% Zhang et al., J Zhejiang Univ Eng Sci54(2020)118, Eqs6,9 and Tables5,8,9.
% Uses original coefficients without fitting: thermal-oxidative aging only.
% Opening range is restricted to 7..10 mm, the reported waterproof tests.
% Ratios are relative to each opening's unaged gasket. Neither an absolute
% leakage rate nor a freeze-thaw-to-aging equivalence is implied.
arguments
    ageYears double {mustBeReal,mustBeFinite,mustBeNonnegative}
    openingMillimetres double {mustBeReal,mustBeFinite}
end
assert(all(ageYears(:) <= 100), 'GasketReference:Age', ...
    'Published prediction range is at most 100 years; no longer extrapolation.');
assert(all(openingMillimetres(:) >= 7 & openingMillimetres(:) <= 10), ...
    'GasketReference:Opening','Reported waterproof tests use openings of 7, 9 and 10 mm.');
timePower = ageYears.^0.41601;
materialRatio = exp(-0.14859*timePower);
openingExponent = -0.136*openingMillimetres + 2.377;
waterproofRatio = exp(-0.14859*openingExponent.*timePower);
result = struct('materialRetention',materialRatio, ...
    'waterproofPressureRetention',waterproofRatio,'ageYears',ageYears, ...
    'opening_mm',openingMillimetres,'serviceTemperature_K',293, ...
    'doi','10.3785/j.issn.1008-973X.2020.01.014', ...
    'source','https://www.zjujournals.com/eng/article/2020/1008-973X/202001014.shtml', ...
    'scope','Published empirical formula reproduction; accelerated thermal oxidation; no freeze-thaw calibration or independent gasket leakage validation.');
end
