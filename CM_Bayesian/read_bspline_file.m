function [data,shift,rho_factor] = read_bspline_file(filename)
%read_bspline_file Reads B-Spline data exported from python
%   "data" contains columns [knots, Fs] of the bspline
%   "rho_factor" is the rho_factor for the parameterisation
lines = readlines(filename);
for ln = 1:length(lines)
    line = convertStringsToChars(lines(ln));
    if ln == 4
        shift = sscanf(line(8:end),'%e %e');
    elseif ln == 5
        rho_factor = sscanf(line(13:end),'%e');
    elseif ln ==6
        n_knots = sscanf(line(13:end),'%d');
        data = zeros(n_knots,2);
    elseif (ln > 6) && (ln <= n_knots+6)
        data(ln-6,:) = sscanf(line,'%f %e');
    end
end
end