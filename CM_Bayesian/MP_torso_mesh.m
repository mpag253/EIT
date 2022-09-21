

clc, clear all, close

% Set shrink factor (closer to 1 is tighter)
% generally want as high as possible without obvious overshrink

filename = './geom/ST4_A'
shrink_factor = 0.7;

disp([filename,'.stl'])
geom = stlread([filename,'.stl']);

geom_nodes = geom.Points;
x = geom_nodes(:, 1);
y = geom_nodes(:, 2);
z = geom_nodes(:, 3);

boundary_i = boundary(x, y, shrink_factor);
x_boundary = x(boundary_i); 
y_boundary = y(boundary_i);
% [x_boundary, y_boundary]

plot(x, y, '.k', 'MarkerSize', 6), hold on
plot(x_boundary, y_boundary, 'r')
plot(x_boundary, y_boundary, 'or', 'MarkerSize', 6)
xlabel('x'), ylabel('y'), axis equal
title({strcat(filename,'.stl boundary'), ...
       strcat('(shrink factor = ', string(shrink_factor), ')')}, ...
    'Interpreter', 'none', 'FontSize', 14)  
set(gcf, 'units', 'normalized')
set(gcf, 'outerposition', [0 0 1 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate mesh
% mesh_shape = @(p)dcircle(p,0,0,1);                                           
% pfix = [];
% figure(101)
% [Nodes,~] = distmesh2d(mesh_shape,@huniform,h,[-1,-1;1,1],pfix);              % THIS PLOTS
Nodes = [x, y];
dt = delaunayTriangulation(Nodes);
Tris = dt.ConnectivityList;

% All nodes
x = Nodes(:,1); y = Nodes(:,2);                   % cartesian coords of nodes
theta = acos(x);                                % polar angles of nodes
nn = length(Nodes);                               % # nodes

% All boundary nodes
full_bdy = x.^2+y.^2>0.999;                       % boundary nodes            (CHANGE CONDITION FOR UNIVERSAL!!!)
bdy_indx = find(full_bdy);                        % indices of boundary nodes
nn_b = sum(full_bdy);                             % # nodes on boundary
[~,ss_b] = sort(atan2(y(full_bdy),x(full_bdy)));  % indices when theta-sorted

% Robin boundary condition
bdy_rob = full_bdy & abs(x)<=.2 & y<0;            % selection of robin BC
bdy_rob_indx = find(bdy_rob);                     % indices of robin BC
nn_rob = sum(bdy_rob);                            % # nodes with robin BC
[~,ss_rob] = sort(atan2(y(bdy_rob),x(bdy_rob)));  % indices when theta-sorted

% Neumann boundary condition
bdy_neu = full_bdy & ~bdy_rob;                    % neumann BC if not robin
bdy_neu_indx = find(bdy_neu);                     % indices of neumann BC
nn_neu = sum(bdy_neu);                            % # nodes with neumann BC
[~,ss_neu] = sort(atan2(y(bdy_neu),x(bdy_neu)));  % indices when theta-sorted