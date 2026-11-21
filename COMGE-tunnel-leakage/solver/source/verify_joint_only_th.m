function checks=verify_joint_only_th()
%VERIFY_JOINT_ONLY_TH Current v16 small-mesh monolithic verification.
% Exercises freezing, water/heat balance, an actually open defect, mechanical
% feedback and the complete uncompressed-vs-coloured tangent. This is a
% numerical verification, not a calibrated annual or 20-year prediction.
checks=TH_2D_segment_jointonly('selftest');
out=checks.coupledOutput;
assert(strcmp(out.version,'THM_v16_20260930_segment_contact'));
assert(out.parameters.liningK==0&&out.leakage.bodyVolume_m3_per_m==0);
assert(all(out.leakage.volumeByJoint_m3_per_m([1:6,8])==0));
checks.impermeableConcrete=true;
checks.sealedJointsStayImpermeable=true;
disp(rmfield(checks,{'coupledOutput','uncoupledOutput'}));
end
