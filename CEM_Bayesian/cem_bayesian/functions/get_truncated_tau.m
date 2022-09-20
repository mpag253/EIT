function [tau_trunc] = get_truncated_tau(tau)
    tau_trunc = tau;
    tau_trunc(tau_trunc>1) = 1;
    tau_trunc(tau_trunc<0) = 0;
end