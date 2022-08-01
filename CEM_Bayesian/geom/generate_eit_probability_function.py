# import openpyxl
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.path as pth
# from matplotlib import cm
# import Slice2D.element_plane_intersection_curve as elemcurve
# import Slice2D.assembly_plane_intersection_curve as assmcurve
# import Slice2D.read_in_data as read_in_data
import BSplineParameterisation.bspline_parameterisation as bsparam
# import scipy.interpolate as ip
# import time
# from sys import exit
# import copy
import pickle
# from scipy import stats
# from scipy import sparse
# import scipy.io
from matplotlib import path as pth
# import parameterise_torso as parat
import os

# INPUTS ###############################################################################################################

pca_id = 'LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A'
mesh_id = 'trimesh_nodes_0010'

version = 'sample_mean'
# version = 'predicted'

# Define the subject ID for the Leave-Out-Out
subject_ids = {'A':'H5977', 'B':'AGING043', 'C':'H7395', 'D':'AGING014', 'E':'AGING053'}
subject_id = subject_ids[pca_id.split('-')[-1]]

# Conditional covariance
if version == 'sample_mean':
    Gamma_b = np.load('pca_'+pca_id+'/sample_mean/sample_bsplines_condcov.npy')
    Gamma_b += 1e-6*np.eye(len(Gamma_b))
elif version == 'predicted':
    Gamma_bgivenp = np.load('pca_'+pca_id+'/predicted_'+subject_id+'/predicted_bsplines_condcov.npy')
    Gamma_bgivenp += 1e-6*np.eye(len(Gamma_bgivenp))

# Parameterised shapes
if version == 'sample_mean':
    bspline_data = np.load('pca_'+pca_id+'/sample_mean/sample_bsplines_mean.npy')
elif version == 'predicted':
    bspline_data = np.load('pca_' + pca_id + '/predicted_'+subject_id+'/predicted_bsplines_mean.npy')
[Is_body, Js_body, Fs_mean, Ks_mean] = bspline_data.T
[lung_index_r, lung_index_l] = [0, 1]
lung_indices = [lung_index_r, lung_index_l]
n_evals = 1000
eval_phis = [np.linspace(0, 1, n_evals),
             np.linspace(0, 1, n_evals),
             np.linspace(0, 2*np.pi, n_evals)]

# Bayesian
n_samples = 10000

# Mesh
trimesh_file = 'pca_'+pca_id+'/sample_mean/trimesh/'+mesh_id+'.csv'
trimesh_nodes = np.genfromtxt(trimesh_file, delimiter=',')
trimesh_bspline_data = np.load('pca_'+pca_id+'/sample_mean/sample_bsplines_mean.npy')
[Is_trimesh, Js_trimesh, Fs_trimesh, Ks_trimesh] = bspline_data.T
Fs_trimesh = Fs_trimesh[Is_trimesh==2]
Ks_trimesh = Ks_trimesh[Is_trimesh==2]
n_nodes = len(trimesh_nodes)

# Setup
bodies = np.unique(Is_body).astype(int)
body_dims = [np.unique(Js_body[Is_body==i]).astype(int) for i in bodies]
n_body_dims = np.sum([len(i) for i in body_dims])
show_plots = False


def pickle_dump(filename, data):
    # Save data to file (pickle)
    pkl_file = open(filename + '.pkl', "wb")
    pickle.dump(data, pkl_file)
    pkl_file.close()
    print("\nSaved data to pickle: ", filename + '.pkl')


def test_complete_intcurves(int_curves):
    """Test if the intersection curve is complete"""

    accepted = True
    for i in range(len(int_curves)):
        int_curve = int_curves[i]
        if not np.all(np.isclose(int_curve[0, 2:], int_curve[-1, 2:], atol=1e-01)):
            accepted = False

    return accepted


def distort_mesh(Fs_sbj, Ks_sbj, sbj_mesh_centre=None):
    """distort the mesh nodes based on fitted and imported bsplines"""

    # # fit bspline to subject torso
    # body_coords = int_curves[2][:, 2:]
    # if sbj_mesh_centre is not None:
    #     body_coords[:, :2] -= sbj_mesh_centre
    # # print(body_coords)
    # plt.close('all')
    # knots_sbj = np.linspace(0, 2*np.pi, 100)[1:]
    # [Fs_sbj, _] = bsparam.get_parameterised_curve(body_coords, knots_sbj, convert_polar=True, convert_rho=False,
    #                                               discont_elems=None, plot_bsplines=False, plot_bases=False)

    # angular position of each node
    node_theta = np.arctan2(trimesh_nodes[:, 1], trimesh_nodes[:, 0]) + np.pi
    # radius of each node
    node_base_r = np.sqrt(np.square(trimesh_nodes[:, 1]) + np.square(trimesh_nodes[:, 0]))
    # base torso radius at each node angular position
    node_base_ro = bsparam.evaluate_mesh_bspline(Fs_trimesh, Ks_trimesh, node_theta)
    # sampled torso radius at each node angular position
    # print(Fs_sbj)
    # print(Ks_sbj)
    node_sbj_ro = bsparam.evaluate_mesh_bspline(Fs_sbj, Ks_sbj, node_theta)

    # new node radius
    node_sbj_r = np.multiply(node_base_r, np.divide(node_sbj_ro, node_base_ro))
    # convert back to cartesian
    # node_base_x = -np.multiply(node_base_r, np.cos(node_theta))
    # node_base_y = -np.multiply(node_base_r, np.sin(node_theta))
    # node_base_xo = -np.multiply(node_base_ro, np.cos(node_theta))
    # node_base_yo = -np.multiply(node_base_ro, np.sin(node_theta))
    # node_sbj_xo = -np.multiply(node_sbj_ro, np.cos(node_theta))
    # node_sbj_yo = -np.multiply(node_sbj_ro, np.sin(node_theta))
    node_sbj_x = -np.multiply(node_sbj_r, np.cos(node_theta))
    node_sbj_y = -np.multiply(node_sbj_r, np.sin(node_theta))
    #
    # # plot check
    # fig, (ax1, ax2) = plt.subplots(1, 2)
    # ax1.plot(body_coords[:, 0], body_coords[:, 1], color="blue")  # blue
    # ax1.scatter(node_base_x, node_base_y, color="black", marker='x')         # blue
    # ax1.scatter(node_base_xo, node_base_yo, color="green", marker='.')       # ? FAIL
    # ax2.scatter(node_sbj_xo, node_sbj_yo, color="red", marker='.')          # ? FAIL
    # ax2.scatter(node_sbj_x, node_sbj_y, color="black", marker='x')
    # plt.show()
    #
    subject_mesh_nodes = np.hstack((node_sbj_x[:, None], node_sbj_y[:, None]))
    # print(subject_mesh_nodes)

    # Plot the morph
    ord = np.argsort(node_theta)
    xo_base = np.multiply(node_base_ro[ord], -np.cos(node_theta[ord]))
    yo_base = np.multiply(node_base_ro[ord], -np.sin(node_theta[ord]))
    xo_sbj = np.multiply(node_sbj_ro[ord], -np.cos(node_theta[ord]))
    yo_sbj = np.multiply(node_sbj_ro[ord], -np.sin(node_theta[ord]))

    if show_plots:
        plt.figure(figsize=(1000/96, 800/96), dpi=96)
        plt.scatter(trimesh_nodes[:, 0], trimesh_nodes[:, 1], color='blue', s=3)
        plt.plot(xo_base, yo_base, color='blue')
        plt.xlim([-175, 175])
        plt.ylim([-125, 125])
        plt.gca().set_aspect('equal', adjustable='box')

        plt.figure(figsize=(1000/96, 800/96), dpi=96)
        plt.scatter(trimesh_nodes[:, 0], trimesh_nodes[:, 1], color='blue', s=3)
        plt.scatter(node_sbj_x, node_sbj_y, color='orange', s=3)
        plt.plot(xo_base, yo_base, color='blue')
        plt.plot(xo_sbj, yo_sbj, color='orange')
        plt.quiver(trimesh_nodes[:, 0], trimesh_nodes[:, 1],
                   node_sbj_x - trimesh_nodes[:, 0], node_sbj_y - trimesh_nodes[:, 1],
                   angles='xy', scale_units='xy', scale=1)
        plt.xlim([-175, 175])
        plt.ylim([-125, 125])
        plt.gca().set_aspect('equal', adjustable='box')

        plt.figure(figsize=(1000/96, 800/96), dpi=96)
        plt.scatter(node_sbj_x, node_sbj_y, color='orange', s=3)
        plt.plot(xo_sbj, yo_sbj, color='orange')
        plt.xlim([-175, 175])
        plt.ylim([-125, 125])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        plt.figure(figsize=(1000/96, 800/96), dpi=96)
        plt.scatter(node_sbj_x, node_sbj_y, color='orange', s=3)
        plt.plot(xo_sbj, yo_sbj, color='orange')
        plt.xlim([-175, 175])
        plt.ylim([-125, 125])
        plt.gca().set_aspect('equal', adjustable='box')

    return subject_mesh_nodes


def get_nodes_in_lung(lung_coordinates, subject_mesh_nodes, sbj_mesh_centre=None):
    """Determine if each node is within the specified lung"""

    if sbj_mesh_centre is not None:
        lung_coordinates -= sbj_mesh_centre
    lung_path = pth.Path(lung_coordinates)
    lung_nodes = lung_path.contains_points(subject_mesh_nodes)

    return lung_nodes.astype(int)


def generate_sampled_lung_nodes():

    #print('\n')
    progress_levels = np.arange(10)/10
    progress_level_indices = (n_samples*progress_levels).astype(int)
    samp_lung_nodes_r = np.empty([n_nodes, n_samples])
    samp_lung_nodes_l = np.empty([n_nodes, n_samples])

    # Iterate until accepted samples = n_samples
    for samp in range(n_samples):

        # Generate random sample from normal distribution
        rand_samp = np.random.normal(0, 1, np.shape(Fs_mean))
        # Cholesky factorisation of covariance matrix
        if version == 'sample_mean':
            nL = np.linalg.cholesky(Gamma_b)
        elif version == 'predicted':
            nL = np.linalg.cholesky(Gamma_bgivenp)
        # Randomly sampled bspline data
        bspline_samp = np.matmul(nL, rand_samp)
        Fs_samp = Fs_mean + bspline_samp

        # Evaluate the sampled parameterised shape
        eval_shapes = [[] for _ in bodies]
        for body in bodies:
            eval_shape_body = np.empty([0, len(eval_phis[body])])
            for dim in body_dims[body]:
                body_dim_indices = np.multiply((Is_body == body), (Js_body == dim))
                Fs_body_dim = Fs_samp[body_dim_indices]
                Ks_body_dim = Ks_mean[body_dim_indices]
                eval_shape = bsparam.evaluate_mesh_bspline(Fs_body_dim, Ks_body_dim, eval_phis[body])
                # print(body, dim)
                # print(np.shape(eval_shape))
                eval_shape_body = np.vstack([eval_shape_body, eval_shape])
            # print(np.shape(eval_shape_body))
            eval_shapes[body] = eval_shape_body
            # print(eval_shapes[body])

        # Generate subject specific mesh
        torso_indices = Is_body==2
        Fs_mode_torso = Fs_samp[torso_indices]
        Ks_mean_torso = Ks_mean[torso_indices]
        subject_mesh_nodes = distort_mesh(Fs_mode_torso, Ks_mean_torso)

        # Determine if each node is within the lungs
        lung_r_coords = eval_shapes[0].T
        lung_l_coords = eval_shapes[1].T
        # print(lung_r_coords)
        # print(lung_l_coords)
        samp_lung_nodes_r[:, samp] = get_nodes_in_lung(lung_r_coords, subject_mesh_nodes)
        samp_lung_nodes_l[:, samp] = get_nodes_in_lung(lung_l_coords, subject_mesh_nodes)

        if show_plots:
            plt.plot(lung_r_coords[:, 0], lung_r_coords[:, 1], color='black')
            plt.plot(lung_l_coords[:, 0], lung_l_coords[:, 1], color='black')
            plot_l = subject_mesh_nodes[samp_lung_nodes_r[:, samp] == 1, :]
            plot_r = subject_mesh_nodes[samp_lung_nodes_l[:, samp] == 1, :]
            plt.scatter(plot_l[:, 0], plot_l[:, 1], color='green', s=8)
            plt.scatter(plot_r[:, 0], plot_r[:, 1], color='green', s=8)
            plt.show()
            plt.close('all')

        if samp in progress_level_indices:
            progress = 100*progress_levels[np.where(progress_level_indices==samp)][0]
            print("Generating {:d} lung node samples... {:.0f}%".format(n_samples, progress), end="\r")

    print("Generating {:d} lung node samples... completed.".format(n_samples), end="\n")
    return samp_lung_nodes_r, samp_lung_nodes_l


def main():
    """..."""

    # Bayesian approach
    samp_lung_nodes_r, samp_lung_nodes_l = generate_sampled_lung_nodes()

    # Node based lung function: mean and covariance (separate lungs):
    uniform_addition = 0.5*np.ones([n_nodes, 1])
    lung_nodes = np.concatenate((samp_lung_nodes_r, samp_lung_nodes_l, uniform_addition), axis=1)
    lung_nodal_mean = np.mean(lung_nodes, axis=1)
    lung_nodal_cov = np.cov(lung_nodes)

    # # To save data...
    save_dir = "Geom/pca_{}/{}/".format(pca_id, mesh_id)
    os.makedirs(save_dir, exist_ok=True)
    np.savetxt(save_dir+"nodal_pf_mean_s-{:d}.csv".format(n_samples), lung_nodal_mean, delimiter=",")
    np.savetxt(save_dir+"nodal_pf_condcov_s-{:d}.csv".format( n_samples), lung_nodal_cov, delimiter=",")

    # Plot results
    fig, ax = plt.subplots(subplot_kw={"projection": "3d", "proj_type": "ortho"})
    ax.plot_trisurf(trimesh_nodes[:, 0], trimesh_nodes[:, 1], lung_nodal_mean,
                    cmap='viridis', linewidth=0, antialiased=False)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d", "proj_type": "ortho"})
    ax.plot_trisurf(trimesh_nodes[:, 0], trimesh_nodes[:, 1], lung_nodal_mean,
                    cmap='viridis', linewidth=0, antialiased=False)
    ax.view_init(azim=90, elev=90)
    plt.show()


if __name__ == "__main__":
    main()