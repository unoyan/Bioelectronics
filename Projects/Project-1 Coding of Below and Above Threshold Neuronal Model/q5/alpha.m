function [ An Am Ah ] = alpha(V,T)
%   Returns rate constant in units  per ms (millisecond)
%   Inputs:  V in [mV] and temperature in [deg C]

global Vr

dV = (V - Vr);

phi = 3^((T-6.3)/10);

%An = phi * (eps + 0.10 - 0.01 .* dV) ./ (eps + exp(1 - 0.1 .* dV) - 1);
An = 0.01 * (10-dV) / ( ( exp((10-dV) / 10) )-1);
%Am = phi * (eps + 2.5 - 0.1 .* dV) ./ (eps + exp(2.5 - 0.1 .* dV) - 1);
Am = 0.1 * (25 - dV) / ((exp((25-dV) / 10)) - 1);
%Ah = phi * 0.07 .* exp(-dV ./ 20);
Ah = 0.07 * exp( -dV / 20 );
end



