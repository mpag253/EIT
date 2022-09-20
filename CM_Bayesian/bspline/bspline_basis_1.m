function [basis1] = bspline_basis_1(phi, knots)
    basis1 = ((phi - knots(2))/(knots(4) - knots(2)))*...
             ((phi - knots(2))/(knots(3) - knots(2)));
end