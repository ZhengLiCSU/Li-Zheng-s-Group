function pressure = iapws_sublimation_reference(T)
%IAPWS_SUBLIMATION_REFERENCE Independent R14-08(2011), equation (6), Pa.
% Source: https://iapws.org/technical-guidance/release/MeltSub.download
% Used for verification only; NOT used to force the new EOS to match.
assert(all(T>=50 & T<=273.16),'IAPWS:Domain','Sublimation correlation domain.');
theta=T/273.16;
a=[-21.2144006,27.3203819,-6.10598130];
b=[0.00333333333,1.20666667,1.70333333];
pressure=611.657*exp(sum(a.*theta(:).^b,2)./theta(:));
end
