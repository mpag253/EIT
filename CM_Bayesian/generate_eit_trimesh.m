clear, clc, close all
addpath('distmesh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  INPUTS

pca_id = 'pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A';
version = 'sample_mean';
% version = 'truth';
% % version = 'predicted';

% Version
if pca_id(end) == 'A', subject_id = 'H5977';
elseif pca_id(end) == 'B', subject_id = 'AGING043';
elseif pca_id(end) == 'C', subject_id = 'H7395';
elseif pca_id(end) == 'D', subject_id = 'AGING014';
elseif pca_id(end) == 'E', subject_id = 'AGING053';
end

if strcmp(version,'sample_mean')
    version_id = 'sample_mean';
    mesh_size = 10;
    path = strcat('geom/',pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat(path,'mesh_seeds_torso_',sprintf('%04d',mesh_size),'.csv');
    mesh_lungs = false;

elseif strcmp(version,'truth')
    version_id = strcat('truth_', subject_id);
    mesh_size = 5;
    path = strcat('geom/',pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat(path,'mesh_seeds_torso_',sprintf('%04d',mesh_size),'.csv');
    pv_file_l = strcat(path,'mesh_seeds_lung_left_',sprintf('%04d',mesh_size),'.csv');
    pv_file_r = strcat(path,'mesh_seeds_lung_right_',sprintf('%04d',mesh_size),'.csv');
    mesh_lungs = true;

elseif strcmp(version,'predicted')
    version_id = strcat('predicted_', subject_id);
    mesh_size = 10;
    path = strcat('geom/',pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat(path,'mesh_seeds_torso_',sprintf('%04d',mesh_size),'.csv');
    mesh_lungs = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pv_t = readmatrix(pv_file_t);
if mesh_lungs == true
    pv_l = readmatrix(pv_file_l);
    pv_r = readmatrix(pv_file_r);
    pv_all = [pv_t;pv_l;pv_r];
else
    pv_all = pv_t;
end

bbox = 1.01*[min(pv_t(:,1)) min(pv_t(:,2)); max(pv_t(:,1)) max(pv_t(:,2))];
% [nodes,tris]=distmesh2d(@dpoly,@huniform,mesh_size,bbox,pv,pv);
[nodes,tris]=distmesh2d(@dpoly,@huniform,mesh_size,bbox,pv_all,pv_t);


tgeom = triangulation(tris,nodes);
% tgeom = delaunayTriangulation(Nodes);
% Tris = tgeom.ConnectivityList;

% boundary_indices = boundary(nodes(:,1), nodes(:,2), shrink_factor);
% boundary_operator(boundary_indices) = true;
% boundary_operator = boundary_operator';


save_file = strrep(strrep(pv_file_t, 'mesh_seeds/mesh_seeds_torso', 'trimesh/trimesh'), '.csv', '.mat');
if ~exist(fileparts(save_file), 'dir'), mkdir(fileparts(save_file)), end
save(save_file,'tgeom')
save_file_nodes = strrep(strrep(save_file, 'trimesh/trimesh', 'trimesh/trimesh_nodes'), '.mat', '.csv');
writematrix(tgeom.Points,save_file_nodes)

disp("Done.")
