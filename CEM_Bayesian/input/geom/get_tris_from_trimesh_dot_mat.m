
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A\truth_H5977\trimesh\trimesh_0060_0015";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-B\truth_AGING043\trimesh\trimesh_0060_0015";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-C\truth_H7395\trimesh\trimesh_0060_0015";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-D\truth_AGING014\trimesh\trimesh_0060_0015";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-E\truth_AGING053\trimesh\trimesh_0060_0015";

% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-E\sample_mean\trimesh\trimesh_0080_0020";

% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A\predicted_H5977\trimesh\trimesh_0080_0020";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-B\predicted_AGING043\trimesh\trimesh_0080_0020";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-C\predicted_H7395\trimesh\trimesh_0080_0020";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-D\predicted_AGING014\trimesh\trimesh_0080_0020";
% path = "pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-E\predicted_AGING053\trimesh\trimesh_0080_0020";


matfile = load(strcat(path, '\trimesh.mat'));


writematrix(matfile.tgeom.ConnectivityList, strcat(path, '\trimesh_tris.csv'))