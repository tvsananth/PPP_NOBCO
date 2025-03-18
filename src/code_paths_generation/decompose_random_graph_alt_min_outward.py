import numpy as np
import metis
from compute_fr2lane import *
from metis_Amat_exits import *

from sim_anneal_paths import *
from compute_A_undirected import *
from compute_sols_undir import *

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
    for j in range(parts.shape[0]):
        ned_parts[parts[j]] += 1
    return ned_parts

def func_parts_adaptive_nparts(data_nodes, data_edges_dir, max_cluster_size=16):

    elist = get_elist(data_nodes, data_edges_dir)
    n_parts = 1
    ned_parts = np.array([1000])
    while (ned_parts.max() > max_cluster_size):
          n_parts += 1
          edgecuts, parts = partition_graph(elist, n_parts)
          ned_parts = count_nodes_per_part(data_edges_dir, parts)
    return n_parts, parts, ned_parts, edgecuts

def do_parts_adaptive(save_dir, max_cluster_size):

    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit = load_data_graph(save_dir)

    n_parts, parts, ned_parts, edgecuts = func_parts_adaptive_nparts(data_nodes, data_edges_dir, max_cluster_size=max_cluster_size)

    return data_nodes, data_edges_dir, nodes_demand, nodes_exit, edgecuts, parts, n_parts

def do_sim_anneal_givb(A_p, b_p, x_p_ini, n_samples, do_ini=False):
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

    x_ids_range = np.arange(0, 2*data_edges_dir.shape[0] + n_nodes)
    x_ids_split = []
    x_ids_split_non = []

    for i_np in range(n_parts):
        x_ids_split_i = []
        for j, d_ei in enumerate(data_edges_dir):
            if ((parts[d_ei[0]] == i_np) & (parts[d_ei[1]] == i_np)):
               x_ids_split_i.append(j)
               x_ids_split_i.append(j + data_edges_dir.shape[0])
        for j_node in range(n_nodes):
            if (parts[j_node] == i_np):
                x_ids_split_i.append(2*data_edges_dir.shape[0] + j_node)
        x_ids_split_i = np.sort(np.array(x_ids_split_i))
        x_ids_split_non_i = np.delete(x_ids_range, x_ids_split_i)
        x_ids_split.append(x_ids_split_i)
        x_ids_split_non.append(x_ids_split_non_i)
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
            print(i_np, j_np, n_parts)

    return x_ids_split, x_ids_split_non


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


def compute_sols_alt_min(A, b, x_ids_split, x_ids_split_non, A_split, A_split_non, n_samples_SA=1, x_ini_size=1000, batch_size=10, n_iter=10000):

    x_ini_1000 = np.int8(np.random.randint(2, size=(x_ini_size, A.shape[1])))
    x_fin_1000 = []

    n_batches = x_ini_size//batch_size
    n_samples = n_samples_SA

    t_only_analog = 0
    for fi_1 in range(n_batches):
        x_ini = x_ini_1000[(fi_1)*batch_size: (fi_1+1)*batch_size, :]
        n_ini = batch_size
        ids_ini = np.arange(0, batch_size)
        x_ini_rand = x_ini.copy()
        x_fin = x_ini.copy()
        i = 0
        while ((i <= n_iter) & (n_ini != 0)):
            for i_p in range(len(x_ids_split)):  
                x_p = x_ini[:, x_ids_split[i_p]]
                x_p_non = x_ini[:, x_ids_split_non[i_p]]
                A_p = A_split[i_p]
                b_p = -(np.matmul(A_split_non[i_p], x_p_non.T) - b)
                x_p_res = []
                for i_ini in range(n_ini):
                    feas_sols_i_ini, t_only_analog_i_ini = do_sim_anneal_givb(A_p, b_p[:, i_ini : (i_ini+1)], x_p[i_ini : (i_ini+1), :], n_samples)
                    x_p_res.append(feas_sols_i_ini[0])
                    t_only_analog += t_only_analog_i_ini
                x_p_res = np.array(x_p_res)
                x_ini[:, x_ids_split[i_p]] = x_p_res[:,:]

            chk_conv = np.matmul(A, x_ini[:,:].T) - b
            chk_conv_bool = ((chk_conv.min(axis=0) == 0) & (chk_conv.max(axis=0) == 0))
            ids_conv = np.where(chk_conv_bool == True)[0]
            x_fin[ids_ini[ids_conv], :] = x_ini[ids_conv, :]
            ids_ini = ids_ini[np.where(chk_conv_bool == False)[0]]
            x_ini = x_ini[np.where(chk_conv_bool == False)[0], :]
            n_ini = ids_ini.shape[0]
            i += 1
        for j1 in range(batch_size):
            x_fin_1000.append(x_fin[j1])
        print(fi_1, n_batches)

    return x_fin_1000, t_only_analog


def do_paths_demand_i(i_demand, save_dir, save_dir_decomp, parts, n_parts, x_ini_size, batch_size, n_iter, n_samples_SA):
    
    t_da_ini = time.time()

    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit = load_data_graph(save_dir)
    n_nodes = data_nodes.shape[0]

    node_demand_i = np.array([nodes_demand[i_demand]])

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir_outward(data_nodes, data_edges_dir, nodes_exit)
    edges, nodes, label_exits = get_nodes_exits(data_nodes, data_edges_dir, nodes_exit)
    b_2fr1, b_2fr1_noedges = compute_b_demand_outward(data_nodes, node_demand_i, nodes_exit)

    A = A_mat_2n_noedges.copy()
    b = b_2fr1_noedges.copy()

    x_ids_split, x_ids_split_non = get_x_splits(data_edges_dir, n_nodes, n_parts, parts)
    A_split, A_split_non = get_A_split(A, x_ids_split, x_ids_split_non)

    x_fin_1000, t_only_analog = compute_sols_alt_min(A, b, x_ids_split, x_ids_split_non, A_split, A_split_non, n_samples_SA=n_samples_SA, x_ini_size=x_ini_size, batch_size=batch_size, n_iter=n_iter)
    x_fin_1000 = np.array(x_fin_1000)
    np.savetxt(save_dir_decomp + '/feas_sols_' + str(i_demand) + '.txt', x_fin_1000)

    chk_conv = np.matmul(A, x_fin_1000[:,:].T) - b
    chk_conv_bool = ((chk_conv.min(axis=0) == 0) & (chk_conv.max(axis=0) == 0))
    ids_conv = np.where(chk_conv_bool == True)[0]
    feas_sols_clean = clean_cycles(x_fin_1000[ids_conv,0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i)
    feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)
    np.savetxt(save_dir_decomp + '/feas_sols_sorted' + str(i_demand) + '.txt', feas_sols_clean_uniq)

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

    data_nodes, data_edges_dir, nodes_demand, nodes_exit, edgecuts, parts, n_parts = do_parts_adaptive(save_dir, max_cluster_size)
    n_nodes = data_nodes.shape[0]

    np.savetxt(save_dir_decomp + '/parts_adaptive_outward.txt', parts)
    parts = np.int32(np.genfromtxt(save_dir_decomp + '/parts_adaptive_outward.txt'))
    n_parts = np.unique(parts).shape[0]

    T_digital = 0
    T_analog = 0    
    t_ini = time.time()
    for i_demand in range(nodes_demand.shape[0]):
        t_only_digital, t_only_analog = do_paths_demand_i(i_demand, save_dir, save_dir_decomp, parts, n_parts, x_ini_size, batch_size, n_iter, n_samples_SA)
        T_digital += t_only_digital
        T_analog += t_only_analog
    t_fin = time.time()
    t_sa = t_fin - t_ini
    t_fin_main = time.time()
    np.savetxt(save_dir_decomp + '/times_sa_alt_min.txt', np.array([t_sa]))
    np.savetxt(save_dir_decomp + '/times_paths_full_run.txt', np.array([t_fin_main - t_ini_main]))
    np.savetxt(save_dir_decomp + '/times_paths_only_digital.txt', np.array([T_digital]))
    np.savetxt(save_dir_decomp + '/times_paths_only_analog.txt', np.array([T_analog]))
    
if __name__ == '__main__':
    main_compute()
