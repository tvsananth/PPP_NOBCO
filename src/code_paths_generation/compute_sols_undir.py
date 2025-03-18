import numpy as np
from sim_anneal_paths import *
from compute_A_undirected import *
from clean_cycles_sols import *
from main_k_shortest_paths import *
from main_k_shortest_paths_exit_i import *
import time
import matplotlib, matplotlib.pyplot

def compute_sols_undir(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time):

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir(data_nodes, data_edges_dir, nodes_exit)

    edges, nodes, label_exits = get_nodes_exits(data_nodes, data_edges_dir, nodes_exit)

    for i_demand in range(nodes_demand.shape[0]):
        t_ini = time.time()
        node_demand_i = np.array([nodes_demand[i_demand]])
        b_2fr1, b_2fr1_noedges = compute_b_demand(data_nodes, node_demand_i, nodes_exit)

        A = A_mat_2n_noedges.copy()
        b = b_2fr1_noedges.copy()

        AA = np.dot(A.T, A)
        h = -2.0*np.dot(b.T, A)
        Q = AA + np.diag(h[0])
        offset = np.dot(b.T, b) + 0.0

        Q2 = np.zeros((A.shape[1], A.shape[1]))
        for i in range(data_edges_dir.shape[0]):
            Q2[i, i + data_edges_dir.shape[0]] = 0.5
            Q2[i + data_edges_dir.shape[0], i] = 0.5
        Q = Q + Q2

        # Define Binary Quadratic Model
        bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)
        simAnnSampler = neal.SimulatedAnnealingSampler()
        sampler = simAnnSampler
        response = sampler.sample(bqm, num_reads=n_samples)
        response = response.aggregate()
        filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]
        feas_sols = response.record.sample[filter_idx]

        feas_sols_clean = clean_cycles(feas_sols, edges, nodes, edges, label_exits, node_demand_i)
        feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)


        feas_times = np.matmul(feas_sols_clean_uniq, tij_edges_undir)
        feas_sols_sorted = feas_sols_clean_uniq[np.argsort(feas_times)[0:100]]
        print(feas_sols_sorted.shape)
        t_fin = time.time()
        sols_SA_time[i_demand] = t_fin - t_ini
        np.savetxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '.txt', feas_sols_sorted)
        #import pdb; pdb.set_trace()


def compute_sols_undir_yens_type_int(save_dir, save_dir_def, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time, yens_type_int=1):
    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir(data_nodes, data_edges_dir, nodes_exit)
    for i_demand in range(nodes_demand.shape[0]):
        tij_edges_undir_copy = tij_edges_undir.copy()
        fi_l = np.genfromtxt(save_dir_def + '/feas_sols_sorted_' + str(i_demand) + '.txt', dtype=np.int32)
        node_demand_i = np.array([nodes_demand[i_demand]])
        b_2fr1, b_2fr1_noedges = compute_b_demand(data_nodes, node_demand_i, nodes_exit)
        #import pdb; pdb.set_trace()
        chk_i = (np.matmul(A_mat_2n_noedges, fi_l.T).T - b_2fr1_noedges.T == 0).all(axis=1)
        feas_sols_sorted = fi_l[chk_i]
        fi_t_i = np.matmul(feas_sols_sorted, tij_edges_undir)
        feas_sols_distinct = []
        for i_j in range(feas_sols_sorted.shape[0]):
            fi_t_i = np.matmul(feas_sols_sorted, tij_edges_undir_copy)
            fi_sol_i = feas_sols_sorted[fi_t_i.argmin()]
            feas_sols_distinct.append(fi_sol_i)
            tij_edges_undir_copy[fi_sol_i == 1] = tij_edges_undir_copy[fi_sol_i == 1]*2.0
        feas_sols_distinct = np.unique(np.array(feas_sols_distinct), axis=0)        
        if (yens_type_int == 1):
           np.savetxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '.txt', feas_sols_distinct)
        if (yens_type_int == 2):
           np.savetxt(save_dir + '/feas_sols_k_sort__' + str(i_demand) + '.txt', feas_sols_sorted)
           np.savetxt(save_dir + '/feas_sols_k_distinct_' + str(i_demand) + '.txt', feas_sols_distinct)


def compute_sols_undir_outward(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time):

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir_outward(data_nodes, data_edges_dir, nodes_exit)

    edges, nodes, label_exits = get_nodes_exits(data_nodes, data_edges_dir, nodes_exit)

    for i_demand in range(nodes_demand.shape[0]):
        t_ini = time.time()
        node_demand_i = np.array([nodes_demand[i_demand]])
        b_2fr1, b_2fr1_noedges = compute_b_demand_outward(data_nodes, node_demand_i, nodes_exit)

        A = A_mat_2n_noedges.copy()
        b = b_2fr1_noedges.copy()

        AA = np.dot(A.T, A)
        h = -2.0*np.dot(b.T, A)
        Q = AA + np.diag(h[0])
        offset = np.dot(b.T, b) + 0.0

        # Define Binary Quadratic Model
        bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)
        simAnnSampler = neal.SimulatedAnnealingSampler()
        sampler = simAnnSampler
        response = sampler.sample(bqm, num_reads=n_samples)
        response = response.aggregate()
        filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]
        feas_sols = response.record.sample[filter_idx]

        feas_sols_clean = clean_cycles(feas_sols[:, 0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i)
        feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)


        feas_times = np.matmul(feas_sols_clean_uniq, tij_edges_undir)
        feas_sols_sorted = feas_sols_clean_uniq[np.argsort(feas_times)[0:100]]
        print(feas_sols_sorted.shape)
        t_fin = time.time()
        sols_SA_time[i_demand] = t_fin - t_ini
        np.savetxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '.txt', feas_sols_sorted)


def compute_sols_yens_k_shortest(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_dir, n_samples, sols_SA_time, yens_type_int=0):
    
    edges, nodes, label_exits = get_nodes_exits(data_nodes, data_edges_dir, nodes_exit)

    for i_demand in range(nodes_demand.shape[0]):
        t_ini = time.time()
        node_demand_i = np.array([nodes_demand[i_demand]])

        if ((yens_type_int == 0) | (yens_type_int == 1)):        
           if (yens_type_int == 0):        
              feas_sols_sorted, a_paths, a_paths_dist = get_feas_sols_yens(nodes, nodes_exit, data_edges_dir, tij_edges_dir, node_demand_i[0], n_samples)

           if (yens_type_int == 1):
              feas_sols_sorted, a_paths, a_paths_dist = get_feas_sols_yens_distinct(nodes, nodes_exit, data_edges_dir, tij_edges_dir, node_demand_i[0], n_samples)
        
           print(feas_sols_sorted.shape)
           t_fin = time.time()
           sols_SA_time[i_demand] = t_fin - t_ini
           np.savetxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '.txt', feas_sols_sorted)

        if (yens_type_int == 2):
           for e_i, exit_node in enumerate(nodes_exit):
               nodes_exit_i = np.array([exit_node])
               feas_sols_sorted1, a_paths, a_paths_dist = get_feas_sols_yens_exit_i(nodes, nodes_exit, nodes_exit_i, data_edges_dir, tij_edges_dir, node_demand_i[0], n_samples)
               np.savetxt(save_dir + '/feas_sols_yens_' + str(i_demand) + '_exit' + str(e_i) + '.txt', feas_sols_sorted1)
               feas_sols_sorted2, a_paths, a_paths_dist = get_feas_sols_yens_distinct_exit_i(nodes, nodes_exit, nodes_exit_i, data_edges_dir, tij_edges_dir, node_demand_i[0], n_samples)
               np.savetxt(save_dir + '/feas_sols_yens_distinct_' + str(i_demand) + '_exit' + str(e_i) + '.txt', feas_sols_sorted2)
               feas_sols_sorted_ei = do_grasp(feas_sols_sorted1, feas_sols_sorted2, tij_edges_dir)
               np.savetxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '_exit' + str(e_i) + '.txt', feas_sols_sorted_ei)

def plot_paths(feas_sols, id_sol, data_nodes, data_edges_dir, nodes_demand, nodes_exit, show_all_edges=False, colc1='r', colc2='r'):

    u_fac = 1.0
    for i in range(data_nodes.shape[0]):
           matplotlib.pyplot.scatter(data_nodes[i,1]*u_fac, data_nodes[i,2]*u_fac, color='b')
           matplotlib.pyplot.annotate(np.int32(data_nodes[i,0]), ([data_nodes[i,1]*u_fac, data_nodes[i,2]*u_fac]))
    if (show_all_edges):
       for i in range(data_edges_dir.shape[0]):
           matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='k', alpha=0.3)
       for i, ei in enumerate(data_edges_dir):
           matplotlib.pyplot.arrow(data_nodes[ei[0],1], data_nodes[ei[0],2], 0.55*(data_nodes[ei[1],1] - data_nodes[ei[0],1]), 0.55*(data_nodes[ei[1],2] - data_nodes[ei[0],2]), shape='full', lw=0, length_includes_head=False, head_width=0.01, color='r')

    for i, ni in enumerate(nodes_exit):
        matplotlib.pyplot.scatter(data_nodes[ni, 1], data_nodes[ni, 2], color='m')

    for i, ni in enumerate(nodes_demand):
        matplotlib.pyplot.scatter(data_nodes[ni, 1], data_nodes[ni, 2], color='g')
    cfl_sol = feas_sols[id_sol].copy()
    cfl_sol[cfl_sol >= 1] = 1
    for i in range(data_edges_dir.shape[0]):
        if (cfl_sol[i] == 1):
           matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color=colc1, ls='-')
        if (cfl_sol[i + data_edges_dir.shape[0]] == 1):
           matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color=colc2, ls='-')
        
    matplotlib.pyplot.show()
