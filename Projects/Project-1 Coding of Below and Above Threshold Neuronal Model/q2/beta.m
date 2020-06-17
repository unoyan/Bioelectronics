function [ Bn Bm Bh ] = beta(V,T)
%   Returns rate constant in units  per ms (millisecond)
%   Inputs:  V in [mV] and temperature in [deg C]

global Vr

dV = (V - Vr);
phi = 3^((T-6.3)/10);

%Bn = phi * 0.125 .* exp(-dV ./ 80);
Bn =  0.125 * exp( -dV / 80 );
%Bm = phi * 4 .* exp(-dV/18);
Bm =  4 * exp( -dV / 18 );
%Bh = phi * 1 ./ (exp(3.0 - 0.1 .* dV) + 1);
Bh = 1 / ( (exp((30 - dV) / 10)) + 1 );
end



