function [basis3] = bspline_basis_3(phi, knots)
    basis3 = ((knots(3) - phi)/(knots(3) - knots(1)))*...
             ((knots(3) - phi)/(knots(3) - knots(2)));
end