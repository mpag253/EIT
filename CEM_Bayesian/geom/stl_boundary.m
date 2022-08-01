
clc, clear all, close

filenames = {'ST4_A', };

% Set shrink factor
% closer to 1 is tighter
% generally want as high as possible without obvious overshrink
shrink_factors = {0.7, }; 

export = false;
% export = true;

for i = 1:length(filenames)

    clf

    filename = filenames{i};
    sf = shrink_factors{i};

    disp([filename,'.stl'])
    geom = stlread([filename,'.stl']);
    geom_nodes = geom.Points;
    
    x = geom_nodes(:, 1);
    y = geom_nodes(:, 2);
    z = geom_nodes(:, 3);
    
    boundary_i = boundary(x, y, sf);
    x_boundary = x(boundary_i); 
    y_boundary = y(boundary_i);
    % [x_boundary, y_boundary]
    
    plot(x, y, '.k', 'MarkerSize', 6), hold on
    plot(x_boundary, y_boundary, 'r')
    plot(x_boundary, y_boundary, 'or', 'MarkerSize', 6)
    xlabel('x'), ylabel('y'), axis equal
    title({strcat(filename,'.stl boundary'), ...
           strcat('(shrink factor = ', string(sf), ')')}, ...
        'Interpreter', 'none', 'FontSize', 14)  
    set(gcf, 'units', 'normalized')
    set(gcf, 'outerposition', [0 0 1 1])

    if export
        output_suffix = ['_boundary_sf',sprintf('%04d',sf*1e3)];
        print([filename, output_suffix ], '-dpng')
        writematrix([x_boundary, y_boundary], [filename,output_suffix,'.csv']) 
    end

end
