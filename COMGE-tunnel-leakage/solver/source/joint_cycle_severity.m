function [severity,dRate,dDose]=joint_cycle_severity(dose,rate,initialSeverity)
%JOINT_CYCLE_SEVERITY Stable exponential irreversible hydraulic state law.
% Shared by evolution, checkpoint validation and observation calibration.
% Broadcasting supports per-joint column settings and path-point matrices.
validateattributes(dose,{'numeric'},{'real','finite','nonnegative'},mfilename,'dose');
validateattributes(rate,{'numeric'},{'real','finite','nonnegative'},mfilename,'rate');
validateattributes(initialSeverity,{'numeric'},{'real','finite','>=',0,'<=',1}, ...
    mfilename,'initialSeverity');
dose=double(dose);rate=double(rate);initialSeverity=double(initialSeverity);
argument=rate.*dose;
severity=initialSeverity+(1-initialSeverity).*(-expm1(-argument));
if nargout>1
    remaining=(1-initialSeverity).*exp(-argument);
    dRate=remaining.*dose;dDose=remaining.*rate;
end
end
