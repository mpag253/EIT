function [tv_mean,tv_stdv] = get_total_variation(nodes,tris,x,y,sig_new,sigma_samples)
    %  first calculate mesh element areas
    tri_areas = zeros([length(tris) 1]);
    for i = 1:length(tri_areas)
        x1 = nodes(tris(i,1),1); y1 = nodes(tris(i,1),2);
        x2 = nodes(tris(i,2),1); y2 = nodes(tris(i,2),2);
        x3 = nodes(tris(i,3),1); y3 = nodes(tris(i,3),2);
        tri_areas(i) = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    end
    %  then gradient and total variance at each element (relative to mean) 
    %  for each posterior sample
    tot_vars = zeros([size(sigma_samples,2) 1]);
    for ii=1:size(sigma_samples,2)
        [dFx,dFy] = trigradient(tris,x,y,sigma_samples(:,ii)-sig_new,'face'); % at faces!
        tot_vars(ii) = sum((abs(dFx)+abs(dFy)).*tri_areas)/sum(tri_areas);
    end
    tv_mean = mean(tot_vars);
    tv_stdv = std(tot_vars);
end