function [out] = get_sigma(tau,deriv)
    if (~exist('deriv', 'var')), deriv = 0; end

%     sig_max = 0.00066;
%     sig_min = 0.0;
    sig_bgd = 0.00048;
    sig_lng = 0.00006;
%     a = (sig_max+sig_min)/2;
%     b = (sig_max-sig_min)/(2*tanh(gam_lim));
%     d = atanh((sig_bgd-a)/b);
%     c = atanh((sig_lng-a)/b) - d;

    if deriv == 0
%         sig = sig_lng*gamma;
%         sig = exp(gamma);
%         sig = sig_lng*exp(log(sig_bgd/sig_lng)*tau);
        sig = sig_bgd*exp(log(sig_lng/sig_bgd)*tau);
%         sig = 0.00030 - 0.00018*tanh(8*gamma-4);
%         sig = a + b*tanh(c*gamma+d);
        out = sig;
    elseif deriv == 1
%         sig_deriv = sig_lng;
%         sig_deriv = exp(gamma);
%         sig_deriv = sig_lng*log(sig_bgd/sig_lng)*exp(log(sig_bgd/sig_lng)*tau);
        sig_deriv = sig_bgd*log(sig_lng/sig_bgd)*exp(log(sig_lng/sig_bgd)*tau);
%         sig_deriv = -(0.00018*8)*(sech(8*gamma-4).^2);
%         sig_deriv = b*c*(sech(c*gamma+d).^2);
        out = sig_deriv;
    end
end