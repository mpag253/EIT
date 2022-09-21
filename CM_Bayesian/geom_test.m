
addpath('distmesh')

pv = readmatrix("torso_mean_poly_0010.csv")
mesh_size = 10;

bbox = 1.01*[min(pv(:,1)) min(pv(:,2)); max(pv(:,1)) max(pv(:,2))];
[p,t]=distmesh2d(@dpoly,@huniform,mesh_size,bbox,pv,pv);