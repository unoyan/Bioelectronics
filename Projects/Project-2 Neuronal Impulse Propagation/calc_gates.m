function [alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h]=  calc_gates(V_m)
    global V_rest
    global T
    dV = V_m-V_rest;
    phi = 3^((T-6.3)/10);
    alpha_n =phi*0.01 * (10-dV) / ( ( exp((10-dV) / 10) )-1);
    beta_n = phi*0.125 * exp(-dV/80);
    alpha_m = phi*0.1 * (25 - dV) / ((exp((25-dV) / 10)) - 1);
    beta_m = phi*4 * exp(-dV/18);
    alpha_h = phi*0.07 * exp( -dV / 20 );
    beta_h = phi*1 / ((exp((30 - dV)/10)) + 1);
end

