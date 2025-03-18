import numpy as np

def compute_A_mat_undir(data_nodes, data_edges_dir, nodes_exit):

    n_nodes = data_nodes.shape[0]
    A_mat_2n = np.zeros((n_nodes+1, 2*data_edges_dir.shape[0]))
    for i in range(data_edges_dir.shape[0]):
        A_mat_2n[data_edges_dir[i,0], i] =  1
        A_mat_2n[data_edges_dir[i,1], i] = -1
        A_mat_2n[data_edges_dir[i,0], i + data_edges_dir.shape[0]] = -1
        A_mat_2n[data_edges_dir[i,1], i + data_edges_dir.shape[0]] =  1
    n_edges = data_edges_dir.shape[0]
    A_mat_2n[n_nodes, 0: n_edges][np.isin(data_edges_dir[:,0], nodes_exit)] = 1
    A_mat_2n[n_nodes, 0: n_edges][np.isin(data_edges_dir[:,1], nodes_exit)] = 1
    A_mat_2n[n_nodes, n_edges: 2*n_edges][np.isin(data_edges_dir[:,0], nodes_exit)] = 1
    A_mat_2n[n_nodes, n_edges: 2*n_edges][np.isin(data_edges_dir[:,1], nodes_exit)] = 1

    A_mat_2n_noedges = np.delete(A_mat_2n, nodes_exit, axis=0)

    return A_mat_2n, A_mat_2n_noedges


def compute_b_demand(data_nodes, nodes_demand, nodes_exit):

    n_nodes = data_nodes.shape[0]
    b_2fr = np.zeros(n_nodes+1)
    b_2fr[nodes_demand] = 1
    b_2fr[nodes_exit] = -1
    b_2fr[-1] = 1
    b_2fr1 = np.zeros((n_nodes+1, 1))
    b_2fr1[:,0] = b_2fr[:]

    b_2fr1_noedges = np.delete(b_2fr1, nodes_exit, axis=0)

    return b_2fr1, b_2fr1_noedges


# Compute A matrix to ensure only one outgoing edge in a path - for simplifying cycle elimination
def compute_A_mat_undir_outward(data_nodes, data_edges_dir, nodes_exit):

    n_nodes = data_nodes.shape[0]
    A_mat_2n = np.zeros((2*n_nodes + 1, 2*data_edges_dir.shape[0] + n_nodes))
    for i in range(data_edges_dir.shape[0]):
        A_mat_2n[data_edges_dir[i,0], i] =  1
        A_mat_2n[data_edges_dir[i,1], i] = -1
        A_mat_2n[data_edges_dir[i,0], i + data_edges_dir.shape[0]] = -1
        A_mat_2n[data_edges_dir[i,1], i + data_edges_dir.shape[0]] =  1

        A_mat_2n[n_nodes + data_edges_dir[i,0], i] =  1
        A_mat_2n[n_nodes + data_edges_dir[i,1], i + data_edges_dir.shape[0]] =  1

    for i1 in range(n_nodes):
        A_mat_2n[n_nodes + i1, 2*data_edges_dir.shape[0] + i1] = 1
 
    n_edges = data_edges_dir.shape[0]
    A_mat_2n[2*n_nodes, 0: n_edges][np.isin(data_edges_dir[:,0], nodes_exit)] = 1
    A_mat_2n[2*n_nodes, 0: n_edges][np.isin(data_edges_dir[:,1], nodes_exit)] = 1
    A_mat_2n[2*n_nodes, n_edges: 2*n_edges][np.isin(data_edges_dir[:,0], nodes_exit)] = 1
    A_mat_2n[2*n_nodes, n_edges: 2*n_edges][np.isin(data_edges_dir[:,1], nodes_exit)] = 1

    nodes_exit1 = np.hstack([nodes_exit, nodes_exit + n_nodes])
    A_mat_2n_noedges = np.delete(A_mat_2n, nodes_exit1, axis=0)

    return A_mat_2n, A_mat_2n_noedges

# Compute b matrix to ensure only one outgoing edge in a path - for simplifying cycle elimination
def compute_b_demand_outward(data_nodes, nodes_demand, nodes_exit):

    n_nodes = data_nodes.shape[0]
    b_2fr = np.zeros(2*n_nodes + 1)
    b_2fr[nodes_demand] = 1
    b_2fr[nodes_exit] = -1
    b_2fr[n_nodes: 2*n_nodes] = 1
    b_2fr[-1] = 1
    b_2fr1 = np.zeros((2*n_nodes+1, 1))
    b_2fr1[:,0] = b_2fr[:]

    nodes_exit1 = np.hstack([nodes_exit, nodes_exit + n_nodes])
    b_2fr1_noedges = np.delete(b_2fr1, nodes_exit1, axis=0)

    return b_2fr1, b_2fr1_noedges

def get_nodes_exits(data_nodes, data_edges_dir, nodes_exit):

    n_nodes = data_nodes.shape[0]
    nodes = np.arange(n_nodes)
    edges = np.zeros((data_edges_dir.shape[0]*2, 2), dtype=np.int32)
    edges[:data_edges_dir.shape[0], 0] = data_edges_dir[:,0]
    edges[:data_edges_dir.shape[0], 1] = data_edges_dir[:,1]
    edges[data_edges_dir.shape[0]: 2*data_edges_dir.shape[0], 0] = data_edges_dir[:,1]
    edges[data_edges_dir.shape[0]: 2*data_edges_dir.shape[0], 1] = data_edges_dir[:,0]
    label_exits = -1*(np.ones(nodes.shape[0]))
    label_exits[nodes_exit] = 1

    return edges, nodes, label_exits
