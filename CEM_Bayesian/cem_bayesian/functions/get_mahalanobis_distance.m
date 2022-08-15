function [m_dist] = get_mahalanobis_distance(x,y,tris,sig_new,sig_t,G_post)
    % weight mahalanobis distance using element areas
    md_wts = zeros(length(x),1);
    elem_areas = zeros(size(tris,1),1);
    for e = 1:size(tris,1)
        elem_areas(e) = polyarea(x(tris(e,:)),y(tris(e,:)));
    end
    elem_areas = elem_areas/sum(elem_areas);
    for n = 1:length(x)
        nod_elems = [find(tris(:,1)==n);find(tris(:,2)==n);find(tris(:,3)==n)];
        md_wts(n) = sum(elem_areas(nod_elems))/3;
    end
    % figure(999), trisurf(tris,x,y,md_wts)
    % size(md_wts)
    % size((sig_t-sig_new)')
    wtd_diff = md_wts.*(sig_t-sig_new);
    % size(wtd_diff)
    m_dist = (wtd_diff'/G_post*wtd_diff)^0.5;
    % m_dist = ((sig_t-sig_new)/(G_post)*(sig_t-sig_new))).^0.5;
    % p-value: d^2 follows the chi-squared distribution with n degrees of freedom, 
    % where n is the number of dimensions of the normal distribution. n = ?
end