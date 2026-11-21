function pressure = threephase_saturation_pressure(T,p,phase)
%THREEPHASE_SATURATION_PRESSURE Chemical-equilibrium root, absolute Pa.
arguments
    T (1,1) double
    p (1,1) struct = joint_threephase_parameters()
    phase (1,1) double {mustBeMember(phase,[1,2])} = 2
end
lp=fzero(@(x)affinity(x,T,p,phase),[log(0.01),log(1e5)]);
pressure=exp(lp);
end
function A=affinity(logPressure,T,p,phase)
q=joint_threephase_gibbs(exp(logPressure)-p.atmosphericPressure_Pa,T,p);
A=q.g(3)-q.g(phase);
end
