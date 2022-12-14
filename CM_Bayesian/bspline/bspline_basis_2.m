function [basis2] = bspline_basis_2(phi, knots)
    basis2 = (((phi - knots(1))/(knots(3) - knots(1)))*...
             ((knots(3) - phi)/(knots(3) - knots(2))) +...
             ((knots(4) - phi)/(knots(4) - knots(2)))*...
             ((phi - knots(2))/(knots(3) - knots(2))));
end