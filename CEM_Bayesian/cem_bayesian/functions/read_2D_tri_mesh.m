function [tgeom,boundary_operator] = read_2D_tri_mesh(filename,shrink_factor,...
                                        centre_data,show_plot)
%read_2D_tri_mesh Reads the mesh nodes and boundary nodes of a .mat file 
% containing a triangulation object
%   ...

% Set shrink factor (closer to 1 is tighter)
% generally want as high as possible without obvious overshrink
% filename = './geom/ST4_A'
% shrink_factor = 0.7;

data = load(filename);
nodes = -data.tgeom.Points;         % negative to rotate mesh
tris = data.tgeom.ConnectivityList;

if centre_data
    nodes(:,1) = nodes(:,1) - (max(nodes(:,1)) + min(nodes(:,1)))/2;
    nodes(:,2) = nodes(:,2) - (max(nodes(:,2)) + min(nodes(:,2)))/2;
end

% boundary_indices = boundary(nodes(:,1), nodes(:,2), shrink_factor);
boundary_elems = freeBoundary(data.tgeom);
boundary_nodes = boundary_elems(:, 1);
boundary_operator(boundary_nodes) = true;
boundary_operator = boundary_operator';

tgeom = triangulation(tris,nodes);

if show_plot
    boundary_nodes = nodes(boundary_indices,:);
    plot(nodes(:,1), nodes(:,2), '.k', 'MarkerSize', 6), hold on
    plot(boundary_nodes(:,1), boundary_nodes(:,2), 'r')
    plot(boundary_nodes(:,1), boundary_nodes(:,2), 'or', 'MarkerSize', 6)
    xlabel('x'), ylabel('y'), axis equal
    title({strcat(filename,'.stl boundary'), ...
           strcat('(shrink factor = ', string(shrink_factor), ')')}, ...
        'Interpreter', 'none', 'FontSize', 14)  
    set(gcf, 'units', 'normalized')
    set(gcf, 'outerposition', [0 0 1 1])
end

fprintf('Read mesh "%s.mat"\n', filename)

end