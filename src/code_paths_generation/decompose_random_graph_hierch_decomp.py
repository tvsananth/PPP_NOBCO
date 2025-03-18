import numpy as np
import metis
from compute_fr2lane import *
from metis_Amat_exits import *

from sim_anneal_paths import *
from compute_A_undirected import *
from compute_sols_undir import *

import time

import os
import argparse
from pathlib import Path


def get_elist(data_nodes, data_edges_dir):

    elist = []
    for i in range(data_nodes.shape[0]):
        elist.append([])
    for i, e_i in enumerate(data_edges_dir):
        elist[e_i[0]].append(e_i[1])
        elist[e_i[1]].append(e_i[0])
    for i in range(data_nodes.shape[0]):
        elist[i] = list(np.sort(np.array(elist[i])))

    return elist


def partition_graph(elist, n_parts):

    if (n_parts != 1):
       (edgecuts, parts) = metis.part_graph(elist, n_parts, recursive=False, objtype='cut', contig=True)
    else:
       parts = np.int32(np.zeros(data_nodes.shape[0]))
       edgecuts = 0
    parts = np.array(parts)

    return edgecuts, parts 


def count_nodes_per_part(data_edges_dir, parts):

    n_parts = np.unique(parts).shape[0]
    ned_parts = np.zeros(n_parts)
    for i, ei in enumerate(data_edges_dir):
        if (parts[ei[0]] == parts[ei[1]]):
           ned_parts[parts[ei[0]]] += 1
    ned_parts = ned_parts*2
    return ned_parts


def func_parts_adaptive_nparts(data_nodes, data_edges_dir, max_cluster_size=16):

    elist = get_elist(data_nodes, data_edges_dir)
    n_parts = 1
    ned_parts = np.array([1000])
    while (ned_parts.max() > max_cluster_size):
          n_parts += 1
          edgecuts, parts = partition_graph(elist, n_parts)
          ned_parts = count_nodes_per_part(data_edges_dir, parts)
    return n_parts, parts, ned_parts


def func_parts(data_nodes, data_edges_dir, n_parts):

    elist = []
    for i in range(data_nodes.shape[0]):
        elist.append([])
    for i, e_i in enumerate(data_edges_dir):
        elist[e_i[0]].append(e_i[1])
        elist[e_i[1]].append(e_i[0])
    for i in range(data_nodes.shape[0]):
        elist[i] = list(np.sort(np.array(elist[i])))
    if (n_parts != 1):
       (edgecuts, parts) = metis.part_graph(elist, n_parts, recursive=False, objtype='cut', contig=True)
    else:
       parts = np.int32(np.zeros(data_nodes.shape[0]))
       edgecuts = 0 
    parts = np.array(parts)

    return edgecuts, parts


def save_parts_adaptive_nparts(parts, save_dir_decomp):

    save_filename_parts = save_dir_decomp + '/parts_full.txt'
    if (not(Path(save_filename_parts).is_file())):
       np.savetxt(save_filename_parts, parts)
    else:
       parts = np.genfromtxt(save_filename_parts)
    return parts   
    

def get_edges_parts(parts, data_edges_dir):

    data_edges_dir_parts = parts[data_edges_dir]
    data_edges_dir_parts = np.sort(data_edges_dir_parts, axis=1)
    data_edges_dir_parts = data_edges_dir_parts[data_edges_dir_parts[:,0] != data_edges_dir_parts[:,1]]
    data_edges_dir_parts = np.unique(data_edges_dir_parts, axis=0)        

    data_nodes_parts = np.zeros((np.unique(parts).shape[0], 3))
    data_nodes_parts[:,0] = np.arange(np.unique(parts).shape[0])

    return data_nodes_parts, data_edges_dir_parts


def do_parts_hierch_level1(parts, data_edges_dir, n_parts1, save_dir_decomp):

    data_nodes_parts, data_edges_dir_parts = get_edges_parts(parts, data_edges_dir)

    save_filename_parts1 = save_dir_decomp + '/parts_hierch_level1.txt'
    if (not(Path(save_filename_parts1).is_file())):
       edgecuts1, parts1 = func_parts(data_nodes_parts, data_edges_dir_parts, n_parts1)
       np.savetxt(save_filename_parts1, parts1)
    else:
       parts1 = np.genfromtxt(save_dir_decomp + '/parts_hierch_level1.txt')
    
    return parts1


def do_sim_anneal_givb(A_p, b_p, x_p_ini, n_samples, t_only_SA, do_ini=False):
    A = A_p.copy()
    b = b_p.copy()
    AA = np.dot(A.T, A)
    h = -2.0*np.dot(b.T, A)
    Q = AA + np.diag(h[0])
    offset = np.dot(b.T, b) + 0.0

    Q2 = np.zeros((A.shape[1], A.shape[1]))
    for i in range(A.shape[1]//2):
        Q2[i, i + A.shape[1]//2] = 0.5
        Q2[i + A.shape[1]//2, i] = 0.5
    Q = Q + Q2
    # Define Binary Quadratic Model
    bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)
    simAnnSampler = neal.SimulatedAnnealingSampler()
    sampler = simAnnSampler
    time_ini_SA = time.time()
    if (do_ini == True):
       response = sampler.sample(bqm, num_reads=n_samples, initial_states=x_p_ini, initial_states_generator="tile")
    if (do_ini == False):
       response = sampler.sample(bqm, num_reads=n_samples*1)
    time_fin_SA = time.time()
    time_SA = time_fin_SA - time_ini_SA
    filter_idx = np.argsort(response.record.energy)
    feas_sols = response.record.sample[filter_idx]
    return feas_sols, time_SA


def get_x_splits(data_edges_dir, n_nodes, n_parts, parts):

    x_ids_range = np.arange(0, 2*data_edges_dir.shape[0])
    x_ids_split = []
    x_ids_split_non = []
    p_ids_split = np.zeros((n_parts, n_parts), dtype=np.int32)

    cc = 0
    for i_np in range(n_parts):
        x_ids_split_i = []
        for j, d_ei in enumerate(data_edges_dir):
            if ((parts[d_ei[0]] == i_np) & (parts[d_ei[1]] == i_np)):
               x_ids_split_i.append(j)
               x_ids_split_i.append(j + data_edges_dir.shape[0])    
        x_ids_split_i = np.sort(np.array(x_ids_split_i))
        print(i_np, n_parts, x_ids_split_i)
        x_ids_split_non_i = np.delete(x_ids_range, x_ids_split_i)
        x_ids_split.append(x_ids_split_i)
        x_ids_split_non.append(x_ids_split_non_i)
        p_ids_split[i_np][i_np] = cc
        cc += 1
        print(i_np, n_parts)

    for i_np in range(n_parts):
        for j_np in range(i_np+1, n_parts):
            x_ids_split_i = []
            for j, d_ei in enumerate(data_edges_dir):
                if ((parts[d_ei[0]] == i_np) & (parts[d_ei[1]] == j_np)):
                   x_ids_split_i.append(j)
                   x_ids_split_i.append(j + data_edges_dir.shape[0])
                if ((parts[d_ei[0]] == j_np) & (parts[d_ei[1]] == i_np)):
                   x_ids_split_i.append(j)
                   x_ids_split_i.append(j + data_edges_dir.shape[0])
            if (len(x_ids_split_i) > 0):
               x_ids_split_i = np.sort(np.array(x_ids_split_i))
               x_ids_split_non_i = np.delete(x_ids_range, x_ids_split_i)
               x_ids_split.append(x_ids_split_i)
               x_ids_split_non.append(x_ids_split_non_i)
               p_ids_split[i_np][j_np] = cc
               cc += 1
            print(i_np, j_np, n_parts)

    return x_ids_split, x_ids_split_non, p_ids_split


def get_A_split(A, x_ids_split, x_ids_split_non):

    A_split = []
    A_split_non = []
    for i1 in range(len(x_ids_split)):       
         x_ids_split_i = x_ids_split[i1]
         x_ids_split_non_i = x_ids_split_non[i1]
         A_split_i = np.zeros((A.shape[0], x_ids_split_i.shape[0]))
         for j_np, x_id_j_np in enumerate(x_ids_split_i):
             A_split_i[:, j_np] = A[:, x_id_j_np]
         A_split_non_i = np.zeros((A.shape[0], x_ids_split_non_i.shape[0]))
         for j_np, x_id_j_np in enumerate(x_ids_split_non_i):
             A_split_non_i[:, j_np] = A[:, x_id_j_np]
         A_split.append(A_split_i)
         A_split_non.append(A_split_non_i)
         print(i1, len(x_ids_split))

    return A_split, A_split_non


def get_splits_all(data_nodes, data_edges_dir, node_demand_i, node_exit_i, n_parts, parts):

    n_nodes = data_nodes.shape[0]

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir(data_nodes, data_edges_dir, node_exit_i)
    edges, nodes, label_exits = get_nodes_exits(data_nodes, data_edges_dir, node_exit_i)
    b_2fr1, b_2fr1_noedges = compute_b_demand(data_nodes, node_demand_i, node_exit_i)

    n_edges = data_edges_dir.shape[0]
    A_mat_2n[n_nodes, 0: n_edges][np.isin(data_edges_dir[:,0], node_exit_i)] = 0
    A_mat_2n[n_nodes, 0: n_edges][np.isin(data_edges_dir[:,1], node_exit_i)] = 1
    A_mat_2n[n_nodes, n_edges: 2*n_edges][np.isin(data_edges_dir[:,0], node_exit_i)] = 1
    A_mat_2n[n_nodes, n_edges: 2*n_edges][np.isin(data_edges_dir[:,1], node_exit_i)] = 0

    A_mat_2n_noedges = np.delete(A_mat_2n, node_exit_i, axis=0)


    A = A_mat_2n_noedges.copy()
    b = b_2fr1_noedges.copy()

    x_ids_split, x_ids_split_non, p_ids_split = get_x_splits(data_edges_dir, n_nodes, n_parts, parts)
    A_split, A_split_non = get_A_split(A, x_ids_split, x_ids_split_non)

    return A_mat_2n, A_mat_2n_noedges, edges, nodes, label_exits, b_2fr1, b_2fr1_noedges, A, b, x_ids_split, x_ids_split_non, A_split, A_split_non, p_ids_split


def get_edges_feas(feas_sol, data_edges_dir, p_ids_split):
    n_edges = data_edges_dir.shape[0]

    edges_non0 = []
    for i in range(n_edges):
        if ((feas_sol[i] == 1) | (feas_sol[i + n_edges] == 1)):
           edges_non0.append([data_edges_dir[i][0], data_edges_dir[i][1]])

    nodes_non0_uniq = np.unique(np.concatenate(edges_non0))
    for i in range(nodes_non0_uniq.shape[0]):
        edges_non0.append([nodes_non0_uniq[i], nodes_non0_uniq[i]])
    edges_non0 = np.array(edges_non0)

    par_ids_split = np.zeros(edges_non0.shape[0], dtype=np.int32)
    for i in range(edges_non0.shape[0]):
        par_ids_split[i] = p_ids_split[edges_non0[i][0]][edges_non0[i][1]]
    par_ids_split = np.unique(par_ids_split)
    return edges_non0, par_ids_split


def get_edges_feas_reorder(feas_sol, data_edges_dir, p_ids_split, node_d_i, node_e_i):
    n_edges = data_edges_dir.shape[0]

    data_edges_undir = np.vstack([data_edges_dir, data_edges_dir])
    data_edges_undir[n_edges: 2*n_edges, 0] = data_edges_dir[0: n_edges, 1]
    data_edges_undir[n_edges: 2*n_edges, 1] = data_edges_dir[0: n_edges, 0]
    edges_non0_undir = data_edges_undir[feas_sol == 1]

    edges_non0_auto = []
    edges_non0_cross = []
    i0_arg = np.where(edges_non0_undir[:,0] == node_d_i)[0][0]    
    while (edges_non0_undir[i0_arg][1] != node_e_i):
          edges_non0_auto.append([edges_non0_undir[i0_arg][0], edges_non0_undir[i0_arg][0]])
          if (edges_non0_undir[i0_arg][0] < edges_non0_undir[i0_arg][1]):
             edges_non0_cross.append([edges_non0_undir[i0_arg][0], edges_non0_undir[i0_arg][1]])
          if (edges_non0_undir[i0_arg][1] < edges_non0_undir[i0_arg][0]):
             edges_non0_cross.append([edges_non0_undir[i0_arg][1], edges_non0_undir[i0_arg][0]])
          i0_arg = np.where(edges_non0_undir[:,0] == edges_non0_undir[i0_arg][1])[0][0]
    edges_non0_auto.append([edges_non0_undir[i0_arg][0], edges_non0_undir[i0_arg][0]])
    if (edges_non0_undir[i0_arg][0] < edges_non0_undir[i0_arg][1]):
       edges_non0_cross.append([edges_non0_undir[i0_arg][0], edges_non0_undir[i0_arg][1]])
    if (edges_non0_undir[i0_arg][1] < edges_non0_undir[i0_arg][0]):
       edges_non0_cross.append([edges_non0_undir[i0_arg][1], edges_non0_undir[i0_arg][0]])
    edges_non0_auto.append([edges_non0_undir[i0_arg][1], edges_non0_undir[i0_arg][1]])
    
    edges_non0 = []
    for i in range(len(edges_non0_cross)):
        edges_non0.append(edges_non0_auto[i])
        edges_non0.append(edges_non0_cross[i])
    edges_non0.append(edges_non0_auto[-1])
    edges_non0 = np.array(edges_non0)

    par_ids_split = np.zeros(edges_non0.shape[0], dtype=np.int32)
    for i in range(edges_non0.shape[0]):
        par_ids_split[i] = p_ids_split[edges_non0[i][0]][edges_non0[i][1]]
    return edges_non0, par_ids_split


def get_edges_feas_reorder_1cluster(data_edges_dir, p_ids_split, node_d_i):
    n_edges = data_edges_dir.shape[0]

    data_edges_undir = np.vstack([data_edges_dir, data_edges_dir])
    data_edges_undir[n_edges: 2*n_edges, 0] = data_edges_dir[0: n_edges, 1]
    data_edges_undir[n_edges: 2*n_edges, 1] = data_edges_dir[0: n_edges, 0]

    edges_non0 = []
    edges_non0.append([node_d_i[0], node_d_i[0]])
    edges_non0 = np.array(edges_non0)

    par_ids_split = np.zeros(edges_non0.shape[0], dtype=np.int32)
    for i in range(edges_non0.shape[0]):
        par_ids_split[i] = p_ids_split[edges_non0[i][0]][edges_non0[i][1]]
    return edges_non0, par_ids_split


def cut_exits(x_ids_split_e, node_exits, i_exit_e, data_edges_dir):
    node_exits_cut_e = np.delete(node_exits, i_exit_e)
    data_vstack = np.vstack([data_edges_dir, data_edges_dir])
    x_cut_e_split = []
    
    for i, x_e in enumerate(x_ids_split_e):
        if (~(np.isin(data_vstack[x_e][0], node_exits_cut_e) | np.isin(data_vstack[x_e][1], node_exits_cut_e))):
           x_cut_e_split.append(x_e)
    x_cut_e_split = np.array(x_cut_e_split)
    x_ids_range = np.arange(0, 2*data_edges_dir.shape[0])
    if (len(x_cut_e_split) != 0):
       x_cut_e_split_non = np.delete(x_ids_range, x_cut_e_split)
    else:
       x_cut_e_split = np.array([], dtype=np.int32)
       x_cut_e_split_non = np.arange(0, 2*data_edges_dir.shape[0])
    return x_cut_e_split, x_cut_e_split_non

    
def compute_sols_alt_min(num_ensemble, A_p, b_p, x_ids_split, x_ids_split_non, A_split, A_split_non, t_only_SA, n_samples_SA=1, x_ini_size=1000, batch_size=10, n_iter=10000, do_SA_direct=False, bool_x_zero=None):
    if (do_SA_direct):
       x_ini_1000 = np.int8(np.random.randint(2, size=(n_samples_SA, A_p.shape[1]))) 
       feas_sols_i_ini, t_only_analog = do_sim_anneal_givb(A_p, b_p, x_ini_1000, n_samples_SA, t_only_SA)
    chk_conv = np.matmul(A_p, feas_sols_i_ini[:,:].T) - b_p
    chk_conv_bool = ((chk_conv.min(axis=0) == 0) & (chk_conv.max(axis=0) == 0))
    ids_conv = np.where(chk_conv_bool == True)[0]
    x_fin_1000 = feas_sols_i_ini[ids_conv]
    is_sol_found = (x_fin_1000.shape[0] > 0)
    return x_fin_1000, is_sol_found, t_only_analog


def do_paths_cluster(node_demand_i_par, node_exit_i_par, save_dir, save_dir_decomp, n_parts, n_parts1, num_ensemble, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct, t_only_SA):

    t_da_ini = time.time()

    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit = load_data_graph(save_dir)
    parts = np.genfromtxt(save_dir_decomp + '/parts_full.txt')
    data_nodes_par, data_edges_dir_par = get_edges_parts(parts, data_edges_dir)
    parts1 = np.genfromtxt(save_dir_decomp + '/parts_hierch_level1.txt')

    data_edges_dir_par = np.int32(data_edges_dir_par)
    n_nodes = data_nodes.shape[0]

    A_mat_2n_par, A_mat_2n_noedges_par, edges_par, nodes_par, label_exits_par, b_2fr1_par, b_2fr1_noedges_par, A_par, b_par, x_ids_split_par, x_ids_split_non_par, A_split_par, A_split_non_par, p_ids_split_par = get_splits_all(data_nodes_par, data_edges_dir_par, node_demand_i_par, node_exit_i_par, n_parts1, parts1)

    x_fin_1000_par, is_sol_found_cluster_par, t_only_analog = compute_sols_alt_min(num_ensemble, A_par, b_par, x_ids_split_par, x_ids_split_non_par, A_split_par, A_split_non_par, t_only_SA, n_samples_SA=n_samples_SA, x_ini_size=x_ini_size, batch_size=batch_size, n_iter=n_iter, do_SA_direct=do_SA_direct)
    x_fin_1000_par = np.array(x_fin_1000_par)

    chk_conv_par = np.matmul(A_par, x_fin_1000_par[:,:].T) - b_par
    chk_conv_bool_par = ((chk_conv_par.min(axis=0) == 0) & (chk_conv_par.max(axis=0) == 0))
    ids_conv_par = np.where(chk_conv_bool_par == True)[0]
    feas_sols_clean_par = clean_cycles(x_fin_1000_par[ids_conv_par,0:(2*data_edges_dir_par.shape[0])], edges_par, nodes_par, edges_par, label_exits_par, node_demand_i_par)
    feas_sols_clean_uniq_par = np.unique(feas_sols_clean_par, axis=0)
    chk_conv_clean_par = np.matmul(A_par, feas_sols_clean_uniq_par[:,:].T) - b_par
    chk_conv_clean_bool_par = ((chk_conv_clean_par.min(axis=0) == 0) & (chk_conv_clean_par.max(axis=0) == 0))

    feas_sols_clean_uniq_par = feas_sols_clean_uniq_par[chk_conv_clean_bool_par[:], :]

    np.savetxt(save_dir_decomp + '/rand_sols_floq_feas_clusters_parallel' + '_demand' + str(node_demand_i_par[0]) + '_exit' + str(node_exit_i_par[0]) + '.txt', feas_sols_clean_uniq_par)
    t_da_fin = time.time()
    t_only_digital = (t_da_fin - t_da_ini) - t_only_analog
    return feas_sols_clean_uniq_par, t_only_digital, t_only_analog


def do_paths_auto(node_demand_i_par, save_dir_decomp, data_edges_dir, data_edges_dir_par, x_ids_split, p_ids_split, nodes_exit, i_demand, i_exit, num_ensemble, A_mat_2n_noedges, A, b, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct):

    x_fin_1000 = []
    for isol_par in range(0, 1):
        edges_non0_isol_par, par_ids_split_isol = get_edges_feas_reorder_1cluster(data_edges_dir_par, p_ids_split, node_demand_i_par)
        x_vars_uniq_isol = []
        x_ids_split_isol = []
        x_ids_split_non_isol = []
        A_split_isol = []
        A_split_non_isol = []
        for ui_isol in range(len(par_ids_split_isol)):
            x_cut_e_split, x_cut_e_split_non = cut_exits(x_ids_split[par_ids_split_isol[ui_isol]], nodes_exit, i_exit, data_edges_dir)
            x_ids_split_isol.append(x_cut_e_split)
            x_ids_split_non_isol.append(x_cut_e_split_non)
            x_vars_uniq_isol.append(x_cut_e_split)
        x_vars_uniq_isol = np.unique(np.concatenate(np.array(x_vars_uniq_isol)))
        bool_x_zero = np.ones(2*data_edges_dir.shape[0], dtype=np.int32)
        bool_x_zero[x_vars_uniq_isol] = 0
        A_split_isol, A_split_non_isol = get_A_split(A_mat_2n_noedges, x_ids_split_isol, x_ids_split_non_isol)
        x_fin_1000_isol_par = compute_sols_alt_min(num_ensemble, A, b, x_ids_split_isol, x_ids_split_non_isol, A_split_isol, A_split_non_isol, n_samples_SA=n_samples_SA, x_ini_size=x_ini_size//4, batch_size=batch_size, n_iter=n_iter, do_SA_direct=do_SA_direct, bool_x_zero=bool_x_zero)
        for jsol_par in range(len(x_fin_1000_isol_par)):
            x_fin_1000.append(x_fin_1000_isol_par[jsol_par])
    x_fin_1000 = np.array(x_fin_1000)

    return x_fin_1000


def do_paths_full_auto_cluster(i_demand, i_exit, node_demand_i_par, node_exit_i_par, save_dir, save_dir_decomp, n_parts, n_parts1, num_ensemble, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct, t_only_SA):
    t_da_ini = time.time()

    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit = load_data_graph(save_dir)
    parts = np.genfromtxt(save_dir_decomp + '/parts_full.txt')
    data_nodes_par, data_edges_dir_par = get_edges_parts(parts, data_edges_dir)
    parts1 = np.genfromtxt(save_dir_decomp + '/parts_hierch_level1.txt')
    data_edges_dir_par = np.int32(data_edges_dir_par)
    n_nodes = data_nodes.shape[0]

    node_demand_i = np.array([nodes_demand[i_demand]])
    node_exit_i = np.array([nodes_exit[i_exit]])

    A_mat_2n, A_mat_2n_noedges, edges, nodes, label_exits, b_2fr1, b_2fr1_noedges, A, b, x_ids_split, x_ids_split_non, A_split, A_split_non, p_ids_split = get_splits_all(data_nodes, data_edges_dir, node_demand_i, nodes_exit, n_parts, parts)

    node_demand_ui_cs = np.array([nodes_demand[i_demand[0]]])
    node_exit_ui_cs = np.array([nodes_exit[i_exit[0]]])

    data_nodes_ui_cs = np.sort(nodes[parts[nodes] == node_demand_i_par])

    fdata_nodes = np.zeros((data_nodes_ui_cs.shape[0], 3))
    fdata_nodes[:, 0] = np.arange(0, data_nodes_ui_cs.shape[0])
    fdata_nodes[:, 1] = data_nodes[data_nodes_ui_cs, 1]
    fdata_nodes[:, 2] = data_nodes[data_nodes_ui_cs, 2]

    fx_ids_split_ui_cs = np.arange(0, data_edges_dir.shape[0])[(parts[data_edges_dir[:, 0]] == node_demand_i_par) & (parts[data_edges_dir[:, 1]] == node_demand_i_par)]

    fdata_edges_dir = []
    fdata_edges_dir_base = []
    for u1, e_u1 in enumerate(fx_ids_split_ui_cs):
        fdata_edges_dir.append([np.where(data_nodes_ui_cs == data_edges_dir[e_u1][0])[0][0], np.where(data_nodes_ui_cs == data_edges_dir[e_u1][1])[0][0]])
        fdata_edges_dir_base.append([data_edges_dir[e_u1][0], data_edges_dir[e_u1][1]])
    fdata_edges_dir = np.array(fdata_edges_dir)
    fdata_edges_dir_base = np.array(fdata_edges_dir_base)
            
    fnode_demand_i = np.where(data_nodes_ui_cs == node_demand_ui_cs)[0]
    fnodes_exit = np.searchsorted(data_nodes_ui_cs, node_exit_ui_cs)

    fn_parts = 1
    fparts = np.zeros(fdata_nodes.shape[0], dtype=np.int32)

    fA_mat_2n, fA_mat_2n_noedges, fedges, fnodes, flabel_exits, fb_2fr1, fb_2fr1_noedges, fA, fb, fx_ids_split, fx_ids_split_non, fA_split, fA_split_non, fp_ids_split = get_splits_all(fdata_nodes, fdata_edges_dir, fnode_demand_i, fnodes_exit, fn_parts, fparts)

    fx_fin_1000_isol_par, is_sol_found_cluster, t_only_analog = compute_sols_alt_min(num_ensemble, fA, fb, fx_ids_split, fx_ids_split_non, fA_split, fA_split_non, t_only_SA, n_samples_SA=n_samples_SA, x_ini_size=x_ini_size//16, batch_size=batch_size, n_iter=n_iter//100, do_SA_direct=do_SA_direct)

    fx_ids_split_ui_cs_undir = np.int32(np.hstack([fx_ids_split_ui_cs, fx_ids_split_ui_cs + data_edges_dir.shape[0]]))
    x_fin_1000 = np.zeros((fx_fin_1000_isol_par.shape[0], 2*data_edges_dir.shape[0]), dtype=np.int32)
    x_fin_1000[:, fx_ids_split_ui_cs_undir] = fx_fin_1000_isol_par[:,:]
    
    np.savetxt(save_dir_decomp + '/rand_sols_floq_feas_full_parallel' + '_demand' + str(i_demand[0]) + '_exit' + str(i_exit[0]) + '.txt', x_fin_1000)

    chk_conv = np.matmul(A, x_fin_1000[:,:].T) - b
    chk_conv_bool = ((chk_conv.min(axis=0) == 0) & (chk_conv.max(axis=0) == 0))
    ids_conv = np.where(chk_conv_bool == True)[0]
    feas_sols_clean = clean_cycles(x_fin_1000[ids_conv,0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i[0])
    feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)
    chk_conv_clean = np.matmul(A, feas_sols_clean_uniq[:,:].T) - b
    chk_conv_clean_bool = ((chk_conv_clean.min(axis=0) == 0) & (chk_conv_clean.max(axis=0) == 0))
    feas_sols_clean_uniq = feas_sols_clean_uniq[chk_conv_clean_bool]
    print('auto', chk_conv_clean_bool)
    np.savetxt(save_dir_decomp + '/rand_sols_floq_feas_full_parallel' + '_demand' + str(i_demand[0]) + '_exit' + str(i_exit[0]) + '.txt', feas_sols_clean_uniq)

    t_da_fin = time.time()
    t_only_digital = (t_da_fin - t_da_ini) - t_only_analog

    return t_only_digital, t_only_analog


def do_paths_full(i_demand, i_exit, node_demand_i_par, node_exit_i_par, save_dir, save_dir_decomp, n_parts, n_parts1, num_ensemble, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct, t_only_SA):

    t_da_ini = time.time()

    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit = load_data_graph(save_dir)
    parts = np.genfromtxt(save_dir_decomp + '/parts_full.txt')
    data_nodes_par, data_edges_dir_par = get_edges_parts(parts, data_edges_dir)
    parts1 = np.genfromtxt(save_dir_decomp + '/parts_hierch_level1.txt')
    data_edges_dir_par = np.int32(data_edges_dir_par)
    n_nodes = data_nodes.shape[0]

    data_edges_undir = np.vstack([data_edges_dir, data_edges_dir])
    data_edges_undir[data_edges_dir.shape[0]:, 0] = data_edges_dir[:, 1]
    data_edges_undir[data_edges_dir.shape[0]:, 1] = data_edges_dir[:, 0]

    node_demand_i = np.array([nodes_demand[i_demand]])
    node_exit_i = np.array([nodes_exit[i_exit]])

    A_mat_2n, A_mat_2n_noedges, edges, nodes, label_exits, b_2fr1, b_2fr1_noedges, A, b, x_ids_split, x_ids_split_non, A_split, A_split_non, p_ids_split = get_splits_all(data_nodes, data_edges_dir, node_demand_i, nodes_exit, n_parts, parts)
    
    feas_sols_clean_uniq_par_1rep = np.genfromtxt(save_dir_decomp + '/rand_sols_floq_feas_clusters_parallel' + '_demand' + str(node_demand_i_par[0]) + '_exit' + str(node_exit_i_par[0]) + '.txt')
    if (len(feas_sols_clean_uniq_par_1rep.shape) == 1):
       feas_sols_clean_uniq_par_1rep = feas_sols_clean_uniq_par_1rep.reshape((1, feas_sols_clean_uniq_par_1rep.shape[0]))
    feas_sum_par_edges = feas_sols_clean_uniq_par_1rep.sum(axis=1)
    feas_arg_par_ids = np.argsort(feas_sum_par_edges)[0:10]
    feas_sols_clean_uniq_par_1rep = feas_sols_clean_uniq_par_1rep[feas_arg_par_ids]
    feas_sols_clean_uniq_par = []
    for i1 in range(feas_sols_clean_uniq_par_1rep.shape[0]):
        for j1 in range(10):
            feas_sols_clean_uniq_par.append(feas_sols_clean_uniq_par_1rep[i1])
    feas_sols_clean_uniq_par = np.array(feas_sols_clean_uniq_par)

    t_only_analog = 0
    x_fin_1000 = []
    for isol_par in range(0, feas_sols_clean_uniq_par.shape[0]):
        feas_isol_par = feas_sols_clean_uniq_par[isol_par]
        edges_non0_isol_par, par_ids_split_isol = get_edges_feas_reorder(feas_isol_par, data_edges_dir_par, p_ids_split, node_demand_i_par, node_exit_i_par)
        x_vars_uniq_isol = []
        x_ids_split_isol = []
        x_ids_split_non_isol = []
        A_split_isol = []
        A_split_non_isol = []
        for ui_isol in range(0, len(par_ids_split_isol)):
            print(ui_isol, len(par_ids_split_isol))
            x_cut_e_split, x_cut_e_split_non = cut_exits(x_ids_split[par_ids_split_isol[ui_isol]], nodes_exit, i_exit, data_edges_dir)
            x_ids_split_isol.append(x_cut_e_split)
            x_ids_split_non_isol.append(x_cut_e_split_non)
            x_vars_uniq_isol.append(x_cut_e_split)
        x_vars_uniq_isol = np.unique(np.concatenate(x_vars_uniq_isol))
        bool_x_zero = np.ones(2*data_edges_dir.shape[0], dtype=np.int32)
        bool_x_zero[x_vars_uniq_isol] = 0

        x_fin_1000_ui = np.zeros(2*data_edges_dir.shape[0], dtype=np.int32)
        node_demand_ui_cs = np.array([nodes_demand[i_demand[0]]])
        ui_coarse = 0
        while (ui_coarse < (len(par_ids_split_isol)//2 + 1)):
            ui_cs = np.array([ui_coarse*2, ui_coarse*2 + 1])
            if (ui_coarse == len(par_ids_split_isol)//2):
               node_exit_ui_cs = np.array([nodes_exit[i_exit[0]]])
            else:
               node_exit_ui_cs = []
               x_ids_exits_ui_cs = x_ids_split_isol[ui_cs[1]]
               x_ids_exits_ui_cs = x_ids_exits_ui_cs[0: len(x_ids_exits_ui_cs)//2]
               for u1, e_u1 in enumerate(x_ids_exits_ui_cs):
                   if (parts[data_edges_dir[e_u1][0]] == edges_non0_isol_par[ui_cs[0]][0]):
                      node_exit_ui_cs.append(data_edges_dir[e_u1][0])
                   if (parts[data_edges_dir[e_u1][1]] == edges_non0_isol_par[ui_cs[0]][0]):
                      node_exit_ui_cs.append(data_edges_dir[e_u1][1])
               node_exit_ui_cs = np.unique(np.array(node_exit_ui_cs))
               node_exit_ui_cs1 = np.delete(node_exit_ui_cs, np.where(node_exit_ui_cs == node_demand_ui_cs)[0])
               if (len(node_exit_ui_cs1) > 0):
                   node_exit_ui_cs = np.delete(node_exit_ui_cs, np.where(node_exit_ui_cs == node_demand_ui_cs)[0])
            x_ids_split_ui_cs = []
            for l1 in range(0, 1):
                x_ids_split_ui_cs.append(x_ids_split_isol[ui_cs[l1]])
            data_nodes_ui_cs = np.sort(nodes[(parts[nodes] == edges_non0_isol_par[ui_cs[0]][0]) | (parts[nodes] == edges_non0_isol_par[ui_cs[0]][1])])
            fdata_nodes = np.zeros((data_nodes_ui_cs.shape[0], 3))
            fdata_nodes[:, 0] = np.arange(0, data_nodes_ui_cs.shape[0])
            fdata_nodes[:, 1] = data_nodes[data_nodes_ui_cs, 1]
            fdata_nodes[:, 2] = data_nodes[data_nodes_ui_cs, 2]
            fdata_edges_dir = []
            fdata_edges_dir_base = []
            for l1 in range(0, 1):
                fx_ids_split_ui_cs = x_ids_split_isol[ui_cs[l1]]
                fx_ids_split_ui_cs = fx_ids_split_ui_cs[0: len(fx_ids_split_ui_cs)//2]
                for u1, e_u1 in enumerate(fx_ids_split_ui_cs):
                    fdata_edges_dir.append([np.where(data_nodes_ui_cs == data_edges_dir[e_u1][0])[0][0], np.where(data_nodes_ui_cs == data_edges_dir[e_u1][1])[0][0]])
                    fdata_edges_dir_base.append([data_edges_dir[e_u1][0], data_edges_dir[e_u1][1]])
            fdata_edges_dir = np.array(fdata_edges_dir)
            fdata_edges_dir_base = np.array(fdata_edges_dir_base)
            fnode_demand_i = np.where(data_nodes_ui_cs == node_demand_ui_cs)[0]
            fnodes_exit = np.searchsorted(data_nodes_ui_cs, node_exit_ui_cs)

            if (len(fx_ids_split_ui_cs) > 0):           

               fx_ids_ui_cs_noexits = fx_ids_split_ui_cs.copy()
               if (len(np.where(np.isin(fdata_edges_dir, fnodes_exit).all(axis=1))[0]) > 0):
                  fx_ids_ui_cs_noexits = np.delete(fx_ids_split_ui_cs, np.where(np.isin(fdata_edges_dir, fnodes_exit).all(axis=1))[0])
                  fdata_edges_dir_base = np.delete(fdata_edges_dir_base, np.where(np.isin(fdata_edges_dir, fnodes_exit).all(axis=1))[0], axis=0)           
                  fdata_edges_dir = np.delete(fdata_edges_dir, np.where(np.isin(fdata_edges_dir, fnodes_exit).all(axis=1))[0], axis=0)

               fdata_edges_undir_base = np.vstack([fdata_edges_dir_base, fdata_edges_dir_base])
               fdata_edges_undir_base[fdata_edges_dir_base.shape[0]:, 0] = fdata_edges_dir_base[:, 1]
               fdata_edges_undir_base[fdata_edges_dir_base.shape[0]:, 1] = fdata_edges_dir_base[:, 0]

               fx_ids_ui_cs_noexits = np.hstack([fx_ids_ui_cs_noexits, fx_ids_ui_cs_noexits + data_edges_dir.shape[0]])
               fn_parts = 1
               fparts = np.zeros(fdata_nodes.shape[0], dtype=np.int32)
               if (fdata_edges_dir.shape[0] != 0):
                  fA_mat_2n, fA_mat_2n_noedges, fedges, fnodes, flabel_exits, fb_2fr1, fb_2fr1_noedges, fA, fb, fx_ids_split, fx_ids_split_non, fA_split, fA_split_non, fp_ids_split = get_splits_all(fdata_nodes, fdata_edges_dir, fnode_demand_i, fnodes_exit, fn_parts, fparts)
                  x_ini_size = 16
                  print(ui_cs, data_nodes_ui_cs, fdata_edges_dir_base, node_demand_ui_cs, node_exit_ui_cs)
                  fx_fin_1000_isol_par, is_sol_found_cluster, t_only_analog_isol_par = compute_sols_alt_min(num_ensemble, fA, fb, fx_ids_split, fx_ids_split_non, fA_split, fA_split_non, t_only_SA, n_samples_SA=n_samples_SA, x_ini_size=x_ini_size//16, batch_size=batch_size, n_iter=n_iter//100, do_SA_direct=do_SA_direct)
                  t_only_analog += t_only_analog_isol_par
               else:
                  is_sol_found_cluster = False
               null_cluster = False
            else:
               is_sol_found_cluster = True
               null_cluster = True

            if (is_sol_found_cluster):
               if (null_cluster == False):
                  fx_fin_1000_isol_par_unique = np.unique(fx_fin_1000_isol_par, axis=0)
                  sols_range_unique = np.arange(0, fx_fin_1000_isol_par_unique.shape[0]) 
                  fx_fin_1000_sol_choose = fx_fin_1000_isol_par_unique[np.random.choice(sols_range_unique)]
                  fx_edges = fedges[fx_fin_1000_sol_choose == 1]
                  fnode_exit_i = fx_edges[:,1][np.where(np.isin(fx_edges[:, 1], fnodes_exit))[0][0]]
                  node_exit_ui_cs_post = data_nodes_ui_cs[fnode_exit_i]
                  if (node_exit_ui_cs_post != node_demand_ui_cs):
                     x_fin_1000_ui[fx_ids_ui_cs_noexits] = fx_fin_1000_sol_choose[:]
                  if (np.isin(node_demand_ui_cs, node_exit_ui_cs).any()):
                     print(node_demand_ui_cs, node_exit_ui_cs_post)
               if (null_cluster == True):
                  node_exit_ui_cs_post = node_demand_ui_cs
               if (ui_coarse != len(par_ids_split_isol)//2):
                  edges_exits_1 = edges[x_ids_split_isol[ui_cs[1]]]
                  if (len(edges_exits_1) > 0):
                     if (np.isin(node_exit_ui_cs_post, edges_exits_1[:,0])):
                        node_demand_ui_cs = np.random.choice(edges_exits_1[:,1][edges_exits_1[:,0] == node_exit_ui_cs_post])
                        print(node_demand_ui_cs)
                        x_fin_1000_ui[x_ids_split_isol[ui_cs[1]][(edges_exits_1[:,0] == node_exit_ui_cs_post) & (edges_exits_1[:,1] == node_demand_ui_cs)]] = 1
                        ui_coarse += 1
                     else:
                        ui_coarse = len(par_ids_split_isol)//2 + 1
                  else:
                     ui_coarse = len(par_ids_split_isol)//2 + 1
               else: 
                  ui_coarse += 1
            else:
                ui_coarse = len(par_ids_split_isol)//2 + 1
        x_fin_1000.append(x_fin_1000_ui)
        x_fin_1000_uif = []
        x_fin_1000_uif.append(x_fin_1000_ui)
        x_fin_1000_uif = np.array(x_fin_1000_uif)
        chk_conv_ui = np.matmul(A, x_fin_1000_uif[:,:].T) - b
        chk_conv_ui_bool = ((chk_conv_ui.min(axis=0) == 0) & (chk_conv_ui.max(axis=0) == 0))
        ids_conv_ui = np.where(chk_conv_ui_bool == True)[0]
        feas_sols_clean_ui = clean_cycles(x_fin_1000_uif[ids_conv_ui, 0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i[0])
        feas_sols_clean_ui_uniq = np.unique(feas_sols_clean_ui, axis=0)
        chk_conv_clean_ui = np.matmul(A, feas_sols_clean_ui_uniq[:,:].T) - b
        chk_conv_clean_ui_bool = ((chk_conv_clean_ui.min(axis=0) == 0) & (chk_conv_clean_ui.max(axis=0) == 0))
        feas_sols_clean_ui_uniq = feas_sols_clean_ui_uniq[chk_conv_clean_ui_bool]
        print(isol_par, chk_conv_ui_bool, chk_conv_clean_ui_bool)
    x_fin_1000 = np.array(x_fin_1000)
    np.savetxt(save_dir_decomp + '/rand_sols_floq_feas_full_parallel' + '_demand' + str(i_demand[0]) + '_exit' + str(i_exit[0]) + '.txt', x_fin_1000)

    chk_conv = np.matmul(A, x_fin_1000[:,:].T) - b
    chk_conv_bool = ((chk_conv.min(axis=0) == 0) & (chk_conv.max(axis=0) == 0))
    ids_conv = np.where(chk_conv_bool == True)[0]
    feas_sols_clean = clean_cycles(x_fin_1000[ids_conv,0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i[0])
    feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)
    chk_conv_clean = np.matmul(A, feas_sols_clean_uniq[:,:].T) - b
    chk_conv_clean_bool = ((chk_conv_clean.min(axis=0) == 0) & (chk_conv_clean.max(axis=0) == 0))
    feas_sols_clean_uniq = feas_sols_clean_uniq[chk_conv_clean_bool]
    np.savetxt(save_dir_decomp + '/rand_sols_floq_feas_full_parallel' + '_demand' + str(i_demand[0]) + '_exit' + str(i_exit[0]) + '.txt', feas_sols_clean_uniq)

    t_da_fin = time.time()
    t_only_digital = (t_da_fin - t_da_ini) - t_only_analog

    return t_only_digital, t_only_analog


def main_compute():

    t_ini_main = time.time()
    parser = argparse.ArgumentParser()

    parser.add_argument('-nodes_file', '--nodes_file_cinput', type=str, help='input nodes file path')
    parser.add_argument('-edges_file', '--edges_file_cinput', type=str, help='input edges file path')
    parser.add_argument('-path_save_data', '--path_save_data_cinput', type=str, help='output directory path')
    parser.add_argument('-decompose_type', '--decompose_type_cinput', type=str, help='Type of decomposition')
    parser.add_argument('-rand_ini_size', '--rand_ini_size_cinput', type=int, default=100, help='number of random seeds')
    parser.add_argument('-batch_size', '--batch_size_cinput', type=int, default=1, help='batch size for processing')
    parser.add_argument('-n_iter', '--n_iter_cinput', type=int, default=10000, help='maximum number of iterations')
    parser.add_argument('-num_samples_SA', '--num_samples_cinput', type=int, default=1, help='number of simulated annealing samples')
    parser.add_argument('-max_cluster_size', '--max_cluster_size_cinput', type=int, help='Maximum number of spins for QUBO')

    args = parser.parse_args()

    nodes_file                 = args.nodes_file_cinput
    edges_file                 = args.edges_file_cinput
    decompose_type             = args.decompose_type_cinput
    save_dir                   = args.path_save_data_cinput + '_' + decompose_type
    x_ini_size                 = args.rand_ini_size_cinput
    batch_size                 = args.batch_size_cinput
    n_iter                     = args.n_iter_cinput
    n_samples_SA               = args.num_samples_cinput
    max_cluster_size           = args.max_cluster_size_cinput

    save_dir_decomp = save_dir

    do_SA_direct = True

    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit = load_data_graph(save_dir)

    n_parts, parts, ned_parts = func_parts_adaptive_nparts(data_nodes, data_edges_dir, max_cluster_size=max_cluster_size)

    parts = save_parts_adaptive_nparts(parts, save_dir_decomp)

    n_parts1 = 1

    parts1 = do_parts_hierch_level1(parts, data_edges_dir, n_parts1, save_dir_decomp)
    time_start = time.perf_counter()

    num_ensemble = 1000                         # number of initial conditions
    n_samples_SA = 1000

    parts = np.genfromtxt(save_dir_decomp + '/parts_full.txt', dtype=np.int32)
    nodes_demand = np.genfromtxt(save_dir + '/nodes_demand.txt', dtype=np.int32)
    nodes_exit = np.genfromtxt(save_dir + '/nodes_exit.txt', dtype=np.int32)

    n_parts = np.unique(parts).shape[0]
     
    t_ini = time.time()

    nodes_demand_par_uniq = np.unique(parts[nodes_demand])
    nodes_exit_par_uniq = np.unique(parts[nodes_exit])
    T_digital = 0
    T_analog = 0
    t_only_SA_clusters = 0

    for i1_d, node_demand_i_par in enumerate(nodes_demand_par_uniq):
        nodes_exit_par_uniq_delete = np.delete(nodes_exit_par_uniq, np.where(nodes_exit_par_uniq == node_demand_i_par)[0])
        for i1_e, node_exit_i_par in enumerate(nodes_exit_par_uniq_delete):
            t_only_SA = 0
            feas_sol_clean_uniq_par_e_d, t_only_digital, t_only_analog = do_paths_cluster(np.array([node_demand_i_par]), np.array([node_exit_i_par]), save_dir, save_dir_decomp, n_parts, n_parts1, num_ensemble, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct, t_only_SA)
            t_only_SA_clusters += t_only_analog
            T_digital += t_only_digital
            T_analog += t_only_analog
            np.savetxt(save_dir_decomp + '/time_onlySA_par_cluster_' + str(node_demand_i_par) + '_' + str(node_exit_i_par) + '.txt', np.array([t_only_SA]))
    t_fin = time.time()

    np.savetxt(save_dir_decomp + '/times_par_cluster.txt', np.array([t_fin - t_ini]))
    np.savetxt(save_dir_decomp + '/times_onlySA_cluster.txt', np.array([t_only_SA_clusters]))

    for i_demand in range(0, nodes_demand.shape[0]):
        t_ini = time.time()
        i_exits = np.arange(0, nodes_exit.shape[0])
        node_demand_i = nodes_demand[i_demand]

        i_demand1 = np.array([i_demand])
        node_demand_i1 = np.array(nodes_demand[i_demand1])
        node_demand_i1_par = np.array(parts[node_demand_i1])

        nodes_exit_par = parts[nodes_exit]
        i_exits_delete = np.delete(i_exits, np.where(nodes_exit_par == parts[node_demand_i])[0])

        t_only_SA = 0
        for i_pool, i_exit in enumerate(i_exits):
            if (node_demand_i1_par == nodes_exit_par[i_exit]):
               t_only_digital, t_only_analog = do_paths_full_auto_cluster(np.array([i_demand]), np.array([i_exit]), node_demand_i1_par, np.array([nodes_exit_par[i_exit]]), save_dir, save_dir_decomp, n_parts, n_parts1, num_ensemble, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct, t_only_SA)
            else:
               t_only_digital, t_only_analog = do_paths_full(np.array([i_demand]), np.array([i_exit]), node_demand_i1_par, np.array([nodes_exit_par[i_exit]]), save_dir, save_dir_decomp, n_parts, n_parts1, num_ensemble, n_samples_SA, x_ini_size, batch_size, n_iter, do_SA_direct, t_only_SA)
            t_only_SA += t_only_analog   
            T_digital += t_only_digital
            T_analog += t_only_analog
        t_fin = time.time()
        np.savetxt(save_dir_decomp + '/times_par_full_demand' + str(i_demand) + '.txt', np.array([t_fin - t_ini]))
        np.savetxt(save_dir_decomp + '/times_onlySA_full_demand' + str(i_demand) + '.txt', np.array([t_only_SA]))

    for i_demand in range(0, nodes_demand.shape[0]):
        fsols = []
        for i_exit in range(0, nodes_exit.shape[0]):
            fsol_i = np.genfromtxt(save_dir_decomp + '/rand_sols_floq_feas_full_parallel' + '_demand' + str(i_demand) + '_exit' + str(i_exit) + '.txt')
            print(i_exit, fsol_i.shape)
            if (len(fsol_i.shape) == 1):
                if (fsol_i.shape[0] > 0):
                   fsol_i = fsol_i.reshape((1, fsol_i.shape[0]))
            if (len(fsol_i.shape) == 2):
               for j_sol in range(fsol_i.shape[0]):
                   fsols.append(fsol_i[j_sol])
        np.savetxt(save_dir_decomp + '/feas_sols_notsorted_' + str(i_demand) + '.txt', np.array(fsols))
        fsols_uniq = np.unique(np.array(fsols), axis=0)
        np.savetxt(save_dir_decomp + '/feas_sols_sorted' + str(i_demand) + '.txt', np.array(fsols_uniq))
    t_fin_main = time.time()
    np.savetxt(save_dir_decomp + '/times_paths_full_run.txt', np.array([t_fin_main - t_ini_main]))
    np.savetxt(save_dir_decomp + '/times_paths_only_digital.txt', np.array([T_digital]))
    np.savetxt(save_dir_decomp + '/times_paths_only_analog.txt', np.array([T_analog]))
    
if __name__ == '__main__':    
    main_compute()
