import numpy as np
from sim_anneal_paths import *
from compute_A_undirected import *
from clean_cycles_sols import *
from main_k_shortest_paths import *
from main_k_shortest_paths_exit_i import *
import time
import matplotlib, matplotlib.pyplot

def get_sampler_dw():
    # Graph corresponding to D-Wave 2000Q
    qpu = DWaveSampler()
    qpu_edges = qpu.edgelist
    qpu_nodes = qpu.nodelist
    if qpu.solver.id == "Advantage_system4.1":
       print(qpu.solver.id)
       X = dnx.pegasus_graph(16, node_list=qpu_nodes, edge_list=qpu_edges)
       dnx.draw_pegasus(X, node_size=1)
       print('Number of qubits=', len(qpu_nodes))
       print('Number of couplers=', len(qpu_edges))
    return qpu

def compute_sols_undir_outward_DW(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time):

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir_outward(data_nodes, data_edges_dir, nodes_exit)

    edges, nodes, label_exits = get_nodes_exits(data_nodes, data_edges_dir, nodes_exit)

    qpu = get_sampler_dw()

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
        bqm_model = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)
        DWavesampler = EmbeddingComposite(qpu)
        DWaveSamples = DWavesampler.sample(bqm=bqm_model, num_reads=n_samples,
                                   return_embedding=True,
                                    annealing_time=1.0
                                   )
        print(round(DWaveSamples.info['embedding_context']['chain_strength'], 3))
        
        response = DWaveSamples.aggregate()
        filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]
        feas_sols = response.record.sample[filter_idx]
        feas_sols_clean = clean_cycles(feas_sols[:, 0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i)
        feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)
        print(feas_sols.shape, feas_sols_clean.shape, feas_sols_clean_uniq.shape)

        solver_greedy = SteepestDescentSolver()
        DWaveSamples_g = solver_greedy.sample(bqm_model, initial_states=DWaveSamples)
        response = DWaveSamples_g.aggregate()
        filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]
        feas_sols = response.record.sample[filter_idx]
        feas_sols_clean = clean_cycles(feas_sols[:, 0:(2*data_edges_dir.shape[0])], edges, nodes, edges, label_exits, node_demand_i)
        feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)
        print(feas_sols.shape, feas_sols_clean.shape, feas_sols_clean_uniq.shape)

        feas_times = np.matmul(feas_sols_clean_uniq, tij_edges_undir)
        feas_sols_sorted = feas_sols_clean_uniq[np.argsort(feas_times)[0:100]]
        print(feas_sols_sorted.shape)
        t_fin = time.time()
        sols_SA_time[i_demand] = t_fin - t_ini
        np.savetxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '.txt', feas_sols_sorted)

