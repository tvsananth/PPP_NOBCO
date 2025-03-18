import numpy as np
import time
from sim_anneal_paths import *

def compute_cmat_undir(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit):

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir_clust(data_nodes, data_edges_dir, nodes_exit)

    edges, nodes, label_exits = get_nodes_exits_clust(data_nodes, data_edges_dir, nodes_exit)

    if (nodes_demand.shape[0] != 0):
       node_demand_i = np.array([nodes_demand[0]])
    else:
       node_demand_i = np.array([])

    for i_demand in range(0, 1):

        t_ini = time.time()
        b_2fr1, b_2fr1_noedges = compute_b_demand_clust(data_nodes, node_demand_i, nodes_exit)

        A = A_mat_2n_noedges.copy()
        b = b_2fr1_noedges.copy()

        AA = np.dot(A.T, A)
        h = -2.0*np.dot(b.T, A)
        Q = AA + np.diag(h[0])
        offset = np.dot(b.T, b) + 0.0

        J = (Q - np.diag(np.diag(Q)))/4.0
        I = np.ones((Q.shape[0], 1))
        H = np.dot(I.T, Q)/2.0

        Q2 = np.zeros((A.shape[1], A.shape[1]))
        for i in range(data_edges_dir.shape[0]):
            Q2[i, i + data_edges_dir.shape[0]] = 0.5
            Q2[i + data_edges_dir.shape[0], i] = 0.5

        J2 = np.zeros((A.shape[1], A.shape[1]))
        H2 = 0.25*(np.ones((1, Q.shape[0])))
        for i in range(data_edges_dir.shape[0]):
            J2[i, i + data_edges_dir.shape[0]] = 0.25*0.5
            J2[i + data_edges_dir.shape[0], i] = 0.25*0.5

        C_mat = np.zeros((J.shape[0] + 1, J.shape[0] + 1))
        C_mat[:-1,:-1] = J[:, :] + J2[:, :]
        C_mat[-1,:-1] = H/2.0 + H2/2.0
        C_mat[:-1, -1] = H/2.0 + H2/2.0


        J = J + J2
        H = H + H2
        min_E = -np.matmul(b.T, b) - 0.25*(np.matmul(I.T, np.matmul(Q, I))) - 0.25*(np.matmul(I.T, np.matmul(np.diag(np.diag(Q)), I))) - 0.25*0.5*(np.dot(I.T, I))

        Q = Q + Q2
    return Q, J, H, C_mat, min_E[0][0], Q2, J2, H2

def compute_A_mat_undir_clust(data_nodes, data_edges_dir, nodes_exit):

    n_nodes = data_nodes.shape[0]
    A_mat_2n = np.zeros((n_nodes, 2*data_edges_dir.shape[0]))
    for i in range(data_edges_dir.shape[0]):
        if (data_edges_dir[i,0] != -1):
           A_mat_2n[data_edges_dir[i,0], i] =  1
           A_mat_2n[data_edges_dir[i,0], i + data_edges_dir.shape[0]] = -1
        if (data_edges_dir[i,1] != -1):
           A_mat_2n[data_edges_dir[i,1], i] = -1
           A_mat_2n[data_edges_dir[i,1], i + data_edges_dir.shape[0]] =  1
    n_edges = data_edges_dir.shape[0]
    if (nodes_exit.shape[0] != 0):
       A_mat_2n_noedges = np.delete(A_mat_2n, nodes_exit, axis=0)
    else:
       A_mat_2n_noedges = A_mat_2n.copy()
    return A_mat_2n, A_mat_2n_noedges


def compute_b_demand_clust(data_nodes, nodes_demand, nodes_exit):

    n_nodes = data_nodes.shape[0]
    b_2fr = np.zeros(n_nodes)
    if (nodes_demand.shape[0] != 0):
       b_2fr[nodes_demand] = 1
    if (nodes_exit.shape[0] != 0):
       b_2fr[nodes_exit] = -1
    b_2fr1 = np.zeros((n_nodes, 1))
    b_2fr1[:,0] = b_2fr[:]
    if (nodes_exit.shape[0] != 0):
       b_2fr1_noedges = np.delete(b_2fr1, nodes_exit, axis=0)
    else:
       b_2fr1_noedges = b_2fr1.copy()
    return b_2fr1, b_2fr1_noedges


def get_nodes_exits_clust(data_nodes, data_edges_dir, nodes_exit):

    n_nodes = data_nodes.shape[0]
    nodes = np.arange(n_nodes)
    edges = np.zeros((data_edges_dir.shape[0]*2, 2), dtype=np.int32)
    edges[:data_edges_dir.shape[0], 0] = data_edges_dir[:,0]
    edges[:data_edges_dir.shape[0], 1] = data_edges_dir[:,1]
    edges[data_edges_dir.shape[0]: 2*data_edges_dir.shape[0], 0] = data_edges_dir[:,1]
    edges[data_edges_dir.shape[0]: 2*data_edges_dir.shape[0], 1] = data_edges_dir[:,0]
    label_exits = -1*(np.ones(nodes.shape[0]))
    if (nodes_exit.shape[0] != 0):
       label_exits[nodes_exit] = 1

    return edges, nodes, label_exits


def compute_sols_undir_clust(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time, clust_id):

    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir_clust(data_nodes, data_edges_dir, nodes_exit)

    edges, nodes, label_exits = get_nodes_exits_clust(data_nodes, data_edges_dir, nodes_exit)

    if (nodes_demand.shape[0] != 0):
       node_demand_i = np.array([nodes_demand[0]])
    else:
       node_demand_i = np.array([])

    for i_demand in range(0, 1):
        t_ini = time.time()
        b_2fr1, b_2fr1_noedges = compute_b_demand_clust(data_nodes, node_demand_i, nodes_exit)

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

        print(feas_sols.shape)
        feas_sols_clean = feas_sols.copy()
        feas_sols_clean_uniq = np.unique(feas_sols_clean, axis=0)


        feas_times = np.matmul(feas_sols_clean_uniq, tij_edges_undir)
        feas_sols_sorted = feas_sols_clean_uniq[np.argsort(feas_times)] #[0:100]]
        print(feas_sols_sorted.shape)

        t_fin = time.time()
        sols_SA_time[i_demand] = t_fin - t_ini
        np.savetxt(save_dir + '/feas_sols_exits_clust' + str(clust_id) + '.txt', feas_sols_sorted)


def plot_paths_clust(feas_sols, id_sol, data_nodes, data_edges_dir, nodes_demand, nodes_exit, data_nodes_parts, data_bound_ids_f, show_all_edges=False):

    matplotlib.pyplot.figure()
    u_fac = 1.0
    for i in range(data_nodes.shape[0]):
           matplotlib.pyplot.scatter(data_nodes[i,1]*u_fac, data_nodes[i,2]*u_fac, color='b')
           matplotlib.pyplot.annotate(np.int32(data_nodes[i,0]), ([data_nodes[i,1]*u_fac, data_nodes[i,2]*u_fac]))
    if (show_all_edges):
       for i in range(data_edges_dir.shape[0]):
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='k')
           if ((data_edges_dir[i,0] == -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.plot([data_nodes_parts[data_bound_ids_f[i],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes_parts[data_bound_ids_f[i],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='b', ls='--')
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] == -1)):
              matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes_parts[data_bound_ids_f[i],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes_parts[data_bound_ids_f[i],2]*u_fac], color='b', ls='--')

       for i, ei in enumerate(data_edges_dir):
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.arrow(data_nodes[ei[0],1], data_nodes[ei[0],2], 0.55*(data_nodes[ei[1],1] - data_nodes[ei[0],1]), 0.55*(data_nodes[ei[1],2] - data_nodes[ei[0],2]), shape='full', lw=0, length_includes_head=False, head_width=0.1, color='r')

    for i, ni in enumerate(nodes_exit):
        matplotlib.pyplot.scatter(data_nodes[ni, 1], data_nodes[ni, 2], color='m')

    for i, ni in enumerate(nodes_demand):
        matplotlib.pyplot.scatter(data_nodes[ni, 1], data_nodes[ni, 2], color='g')
    cfl_sol = feas_sols[id_sol].copy()
    cfl_sol[cfl_sol >= 1] = 1
    for i in range(data_edges_dir.shape[0]):
        if (cfl_sol[i] == 1):
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='r')
           if ((data_edges_dir[i,0] == -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.plot([data_nodes_parts[data_bound_ids_f[i],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes_parts[data_bound_ids_f[i],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='r')
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] == -1)):
              matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes_parts[data_bound_ids_f[i],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes_parts[data_bound_ids_f[i],2]*u_fac], color='r')
        if (cfl_sol[i + data_edges_dir.shape[0]] == 1):
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='r')
           if ((data_edges_dir[i,0] == -1) & (data_edges_dir[i,1] != -1)):
              matplotlib.pyplot.plot([data_nodes_parts[data_bound_ids_f[i],1]*u_fac, data_nodes[data_edges_dir[i,1],1]*u_fac], [data_nodes_parts[data_bound_ids_f[i],2]*u_fac, data_nodes[data_edges_dir[i,1],2]*u_fac], color='r')
           if ((data_edges_dir[i,0] != -1) & (data_edges_dir[i,1] == -1)):
              matplotlib.pyplot.plot([data_nodes[data_edges_dir[i,0],1]*u_fac, data_nodes_parts[data_bound_ids_f[i],1]*u_fac], [data_nodes[data_edges_dir[i,0],2]*u_fac, data_nodes_parts[data_bound_ids_f[i],2]*u_fac], color='r')

    matplotlib.pyplot.show()

def get_a_mat(data_nodes_parts, data_edges_dir_parts, nodes_demand, nodes_exit, parts_2, clust_id):
    data_nodes_f = data_nodes_parts[parts_2 == clust_id]
    data_nodes_ids = np.int32(data_nodes_f[:,0].copy())
    data_nodes_f[:,0] = np.arange(0, data_nodes_f.shape[0])
    labels = {}
    for i in range(data_nodes_ids.shape[0]):
        labels[np.int32(data_nodes_ids[i])] = i
    nodes_exit_ids = nodes_exit[parts_2[nodes_exit] == clust_id]
    nodes_exit_f = []
    for i in range(nodes_exit_ids.shape[0]):
        nodes_exit_f.append(labels[nodes_exit_ids[i]])
    nodes_exit_f = np.int32(np.array(nodes_exit_f))
    nodes_demand_ids = nodes_demand[parts_2[nodes_demand] == clust_id]
    nodes_demand_f = []
    for i in range(nodes_demand_ids.shape[0]):
        nodes_demand_f.append(labels[nodes_demand_ids[i]])
    nodes_demand_f = np.int32(np.array(nodes_demand_f))

    data_edges_dir_f = []
    data_bound_ids_f = []
    data_edges_ids = []
    data_edges_labels = []
    for i, di in enumerate(data_edges_dir_parts):
        if ((parts_2[di[0]] == clust_id) | (parts_2[di[1]] == clust_id)):
           data_edges_ids.append(i)
           a1 = [-1, -1]
           b1 = -1
           if (parts_2[di[0]] == clust_id):
              a1[0] = labels[di[0]]
           else:
              b1 = di[0]
           if (parts_2[di[1]] == clust_id):
              a1[1] = labels[di[1]]
           else:
              b1 = di[1]
           data_edges_dir_f.append(a1)
           data_bound_ids_f.append(b1)
           data_edges_labels.append(i)
    data_edges_dir_f = np.array(np.int32(data_edges_dir_f))
    data_bound_ids_f = np.array(np.int32(data_bound_ids_f))
    data_edges_labels = np.array(np.int32(data_edges_labels))
    return data_nodes_f, data_edges_dir_f, nodes_demand_f, nodes_exit_f, data_bound_ids_f, data_edges_labels

def all_sols(feas_sols1, data_edges_labels, data_edges_dir_parts):

    feas_sols_all = -1*np.int32(np.ones((feas_sols1.shape[0], data_edges_dir_parts.shape[0]*2)))
    feas_sols_all[:, data_edges_labels[:]] = feas_sols1[:, 0:data_edges_labels.shape[0]]
    feas_sols_all[:, data_edges_labels[:] + data_edges_dir_parts.shape[0]] = feas_sols1[:, data_edges_labels.shape[0]:]
    return feas_sols_all

def merge_sols(fsol1, fsol2):
    fsol_12 = []
    for i in range(fsol1.shape[0]):
        fsol1_i = fsol1[i,:].copy()
        for j in range(fsol2.shape[0]):
            fsol2_j = fsol2[j,:].copy()
            fsol12_ij = fsol1_i.copy()
            fsol12_ij[fsol2_j != -1] = fsol2_j[fsol2_j != -1]
            chk_conc = (fsol12_ij[fsol1_i != -1] == fsol1_i[fsol1_i != -1])
            if (chk_conc.all()):
               fsol_12.append(fsol12_ij)
        print (i, fsol1.shape[0], j, fsol2.shape[0])
    fsol_12 = np.array(fsol_12)
    return fsol_12
