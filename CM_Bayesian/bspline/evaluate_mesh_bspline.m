function [rphi] = evaluate_mesh_bspline(Fs, global_knots, phis)
    % reproduce the splines for a given array of phi (phis)

    % Pre-allocate an array for the soloution
    rphi = zeros(length(phis),1);


    % Iterate through each value of phi
    %for i, phi in enumerate(phis):
    for i = 1:length(phis)
        phi = phis(i);

        % Find the local knots that bound the current value of phi
        [local_knots, j] = get_local_knots(phi, global_knots);

        % Evaluate the value of the spline at phi using the basis functions and fitted weights
        % adding periodicity for matlab
        if j==2,        Fs1 = Fs(j-0);  Fs2 = Fs(j-1);      Fs3 = Fs(end);
        elseif j==1,    Fs1 = Fs(j-0);  Fs2 = Fs(end);      Fs3 = Fs(end-1);
        elseif j==0,    Fs1 = Fs(end);  Fs2 = Fs(end-1);    Fs3 = Fs(end-2);  
        else,           Fs1 = Fs(j-0);  Fs2 = Fs(j-1);      Fs3 = Fs(j-2);
        end
        
        rphi(i) = (Fs1*bspline_basis_1(phi, local_knots) +...
                   Fs2*bspline_basis_2(phi, local_knots) +...
                   Fs3*bspline_basis_3(phi, local_knots));
    end
end