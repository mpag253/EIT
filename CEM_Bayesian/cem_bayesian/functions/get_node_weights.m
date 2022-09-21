function [node_wts] = get_node_weights(x,y,tris)

    % get element areas
    elem_areas = zeros(size(tris,1),1);
    for e = 1:size(tris,1)
        elem_areas(e) = polyarea(x(tris(e,:)),y(tris(e,:)));
    end

    % get node areas 
    node_areas = zeros(length(x),1);
    for n = 1:length(x)
        node_elems = [find(tris(:,1)==n);find(tris(:,2)==n);find(tris(:,3)==n)];
        node_areas(n) = sum(elem_areas(node_elems))/3;
    end

    % get node weights
    % relative area of node to average area
    avg_area = sum(elem_areas)/length(x);
    node_wts = node_areas/avg_area;

%     figure(999), trisurf(tris,x,y,node_wts)
%     size(node_wts)
end