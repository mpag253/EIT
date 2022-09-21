function [tgeom_new,var_r_new] = generate_subject_mesh(tgeom,full_bdy,bsp_file_mesh,bsp_file_subj)
%generate_subject_mesh Generates triangle mesh for a subject torso shape
%    Distorts a population mesh (described by a triangulation obect and 
%    boundary b-spline of the radius, r) to a subject mesh (described by a 
%    boundary b-spline of parameter rho)

% [Nodes,Tris,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
% tgeom = triangulation(Tris,Nodes);
% subject_file = "export_torso_param_test.txt";
% mesh_file = "export_torso_param_mesh.txt";

addpath('bspline')

%%%%% UNPACK AND SET UP

Nodes = tgeom.Points;
Tris = tgeom.ConnectivityList;

% All nodes
x = Nodes(:,1); y = Nodes(:,2);               % cartesian coords of nodes
theta = atan2(y,x);                           % polar angles of nodes
r = (x.^2 + y.^2).^.5;
% nn = length(Nodes);                           % # nodes

% All boundary nodes
% bdy_indx = find(full_bdy);                    % indices of boundary nodes
% nn_b = sum(full_bdy);                         % # nodes on boundary
x_bdy = x(full_bdy); y_bdy = y(full_bdy);
theta_bdy = atan2(y_bdy,x_bdy);               % polar angles of boundary nodes
[~,ss_b] = sort(theta_bdy);                   % indices when theta-sorted

%%%%% EVALUATING ORIGINAL AND NEW MESH

% ORIGINAL radius of each boundary node
r_bdy = [x_bdy.^2 + y_bdy.^2].^.5;

% b-spline for the new mesh (rho)
[data,~,rho_factor] = read_bspline_file(bsp_file_subj);

% b-spline for the population mesh (r)
[data_mesh,~,~] = read_bspline_file(bsp_file_mesh);

% rho for the theta of each boundary node
rho_theta_bdy = evaluate_mesh_bspline(data(:,2), data(:,1), theta_bdy+pi);

% NEW radius, r, for the theta of each boundary node
r_theta_bdy = r_bdy + rho_theta_bdy/rho_factor;


%%%%% ADJUSTING MESH

% evaluate rho at boundary for each mesh theta (including internal nodes)
rho_mesh = evaluate_mesh_bspline(data(:,2), data(:,1), theta+pi);

% evaluate old r at boundary for each mesh theta (including internal nodes)
r_mesh_bdy = evaluate_mesh_bspline(data_mesh(:,2), data_mesh(:,1), theta+pi);

% evaluate new r at boundary for each mesh theta (including internal nodes)
r_new_bdy = r_mesh_bdy + rho_mesh/rho_factor;

% evaluate new r for each node
r_new = r.*(r_new_bdy./r_mesh_bdy);

% evaluate new coordinates and mesh
x_new = r_new.*cos(theta);
y_new = r_new.*sin(theta);
tgeom_new = triangulation(Tris,[x_new y_new]);


%%%%% BAYESIAN APPROACH

% (placeholder!) variances at boundary for each node
var_r_new_bdy = 10*ones(length(r_new_bdy),1);
% actually need covariance matrix of b-spline weights
% boundary variances computed using sum formula
% 𝑉𝑎𝑟(𝑎𝑋+𝑏𝑌+𝑐𝑍)=𝑎^2𝑉𝑎𝑟(𝑋)+𝑏^2𝑉𝑎𝑟(𝑌)+𝑐^2𝑉𝑎𝑟(𝑍)+2𝑎𝑏𝐶𝑜𝑣(𝑋,𝑌)+2𝑎𝑐𝐶𝑜𝑣(𝑋,𝑍)+2𝑏𝑐𝐶𝑜𝑣(𝑌,𝑍)

% variances (in radius) for each node
var_r_new = var_r_new_bdy.*(r_new./r);





%%%%% PLOTTING

% % Plot orignal geometry
% figure(1001)
% triplot(tgeom, 'Color', [1. .8 .8]), hold on
% % plot(phis,rphis/rho_factor)
% polar(theta_bdy(ss_b), r_bdy(ss_b),'r')
% polar(theta_bdy(ss_b), r_theta_bdy(ss_b),'b')
% axis equal
% title('Original mesh (blue = subject, red = population)')

% % Plot new geometry
% figure(1002)
% triplot(tgeom_new, 'Color', [.8 .8 1.]), hold on
% % plot(x_new,y_new,'x'), hold on
% % polar(theta,r_new,'.')
% polar(theta_bdy(ss_b), r_bdy(ss_b),'r')
% polar(theta_bdy(ss_b), r_theta_bdy(ss_b),'b')
% axis equal
% title('New (distorted) mesh (blue = subject, red = population)')

% % check true geom
% figure(1003)
% [Nodes_true,Tris_true,full_bdy_true] = read_2D_stl_mesh('./geom/ST4_E',0.7,true,false);
% tgeom_true = triangulation(Tris_true,Nodes_true);
% triplot(tgeom_true, 'Color', [.8 .8 1.]), hold on
% polar(theta_bdy(ss_b), r_bdy(ss_b),'r')
% polar(theta_bdy(ss_b), r_theta_bdy(ss_b),'b')
% axis equal
% title('True mesh (validation)')

end


% function [basis1] = bspline_basis_1(phi, knots)
%     basis1 = ((phi - knots(2))/(knots(4) - knots(2)))*...
%              ((phi - knots(2))/(knots(3) - knots(2)));
% end
% 
% 
% function [basis2] = bspline_basis_2(phi, knots)
%     basis2 = (((phi - knots(1))/(knots(3) - knots(1)))*...
%              ((knots(3) - phi)/(knots(3) - knots(2))) +...
%              ((knots(4) - phi)/(knots(4) - knots(2)))*...
%              ((phi - knots(2))/(knots(3) - knots(2))));
% end
% 
% 
% function [basis3] = bspline_basis_3(phi, knots)
%     basis3 = ((knots(3) - phi)/(knots(3) - knots(1)))*...
%              ((knots(3) - phi)/(knots(3) - knots(2)));
% end
% 
% 
% function [local_knots, j] = get_local_knots(phi, global_knots)
% 
%     % get indices of global_knots points < phi
%     % js = [p for p, q in enumerate(global_knots) if q < phi]
%     lenm = length(global_knots);
%     js = [];
%     for p = 1:lenm, if global_knots(p) < phi, js = [js p]; end, end
% 
%     if isempty(js)
%         j = 0;
%         % local_knots = np.concatenate((global_knots[-2:] - 2*np.pi, global_knots[0:2]))
%         local_knots = [global_knots(end-1:end)'-2*pi global_knots(1:2)'];
% 
%     else
%         j = js(end);
% 
%         if j == 1
%             % local_knots = np.concatenate(([global_knots[-1] - 2*np.pi], global_knots[0:3]))
%             local_knots = [global_knots(end)-2*pi global_knots(1:3)'];
% 
%         elseif j == (lenm)
%             % local_knots = np.concatenate((global_knots[-2:], global_knots[0:2] + 2*np.pi))
%             local_knots = [global_knots(end-1:end)' global_knots(1:2)'+2*pi];
% 
%         elseif j == (lenm-1)
%             % local_knots = np.concatenate((global_knots[-3:], [global_knots[0] + 2*np.pi]))
%             local_knots = [global_knots(end-2:end)' global_knots(1)+2*pi];
% 
%         else
%             % local_knots = global_knots[j - 1:j + 3]
%             local_knots = global_knots(j-1:j+2);
%         end
%     end
% end
% 
% 
% function [rphi] = evaluate_mesh_bspline(Fs, global_knots, phis)
%     % reproduce the splines for a given array of phi (phis)
% 
%     % Pre-allocate an array for the soloution
%     rphi = zeros(length(phis),1);
% 
% 
%     % Iterate through each value of phi
%     %for i, phi in enumerate(phis):
%     for i = 1:length(phis)
%         phi = phis(i);
% 
%         % Find the local knots that bound the current value of phi
%         [local_knots, j] = get_local_knots(phi, global_knots);
% 
%         % Evaluate the value of the spline at phi using the basis functions and fitted weights
%         % adding periodicity for matlab
%         if j==2,        Fs1 = Fs(j-0);  Fs2 = Fs(j-1);      Fs3 = Fs(end);
%         elseif j==1,    Fs1 = Fs(j-0);  Fs2 = Fs(end);      Fs3 = Fs(end-1);
%         elseif j==0,    Fs1 = Fs(end);  Fs2 = Fs(end-1);    Fs3 = Fs(end-2);
%         else,           Fs1 = Fs(j-0);  Fs2 = Fs(j-1);      Fs3 = Fs(j-2);
%         end
%         rphi(i) = (Fs1*bspline_basis_1(phi, local_knots) +...
%                    Fs2*bspline_basis_2(phi, local_knots) +...
%                    Fs3*bspline_basis_3(phi, local_knots));
%     end
% end



% end