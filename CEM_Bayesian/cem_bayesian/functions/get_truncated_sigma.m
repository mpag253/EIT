function [sig_trunc] = get_truncated_sigma(sig)
    sig_trunc = sig;
    sig_trunc(sig_trunc>1) = 1;
    sig_trunc(sig_trunc<0) = 0;
end