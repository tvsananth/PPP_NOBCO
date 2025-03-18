import numpy as np

   
def get_paths(id_d, verts, edts, label_verts, label_exits, label_edges, path_cur_vert, path_cur_edge, paths_de, cc1, paths_cycle):

    if (label_exits[id_d] == 1):
       paths_de.append(path_cur_vert)
       label_edges[path_cur_edge] = 1
       cc1[0] += 1
       return 
    elif ((label_verts[id_d] == 1) | (len(edts[id_d]) == 0)):
       if (len(edts[id_d]) == 0):
          return
       else:
          paths_cycle.append(path_cur_vert)
          return
    else:
        label_verts[id_d] = 1
        for i in range(len(edts[id_d])):
            path_cur_vert_i = path_cur_vert.copy()
            path_cur_vert_i.append(verts[id_d][i])
            path_cur_edge_i = path_cur_edge.copy()
            path_cur_edge_i.append(edts[id_d][i])
            label_verts_i = label_verts.copy()
            get_paths(verts[id_d][i], verts, edts, label_verts_i, label_exits, label_edges, path_cur_vert_i, path_cur_edge_i, paths_de, cc1, paths_cycle)
        return


def get_paths_all(feas_sols_id, edges, nodes, data_edges_dir, label_exits, nodes_demand):

    edges_non0 = data_edges_dir[feas_sols_id == 1]
    edges_non0_ids = np.arange(data_edges_dir.shape[0])[feas_sols_id == 1]

    label_edges = 0*(np.ones(edges.shape[0], dtype=np.int32))

    verts = {}
    edts = {}

    for i in range(edges_non0.shape[0]):
        verts[edges_non0[i,0]] = []
        edts[edges_non0[i,0]] = []
    for i in range(edges_non0.shape[0]):
        verts[edges_non0[i,0]].append(edges_non0[i,1])
        edts[edges_non0[i,0]].append(edges_non0_ids[i])

    paths_de_all = []
    paths_cycle_all = []
    for i in range(nodes_demand.shape[0]):
        id_d = nodes_demand[i]
        label_verts = -1*(np.ones(nodes.shape[0]))
        path_cur_vert = []
        path_cur_edge = []
        cc1 = np.array([0])
        paths_de = []
        paths_cycle = []
        get_paths(id_d, verts, edts, label_verts, label_exits, label_edges, path_cur_vert, path_cur_edge, paths_de, cc1, paths_cycle)
        paths_de_all.append(paths_de)
        paths_cycle_all.append(paths_cycle)

    return paths_de_all, label_edges, paths_cycle_all

def clean_cycles(feas_sols, edges, nodes, data_edges_dir, label_exits, nodes_demand):

    feas_sols_clean = np.zeros((feas_sols.shape[0], feas_sols.shape[1]), dtype=np.int32)

    for id_sol in range(feas_sols.shape[0]):
        feas_sols_id = feas_sols[id_sol]
        paths_de_id, feas_sols_clean_id, paths_cycle_id = get_paths_all(feas_sols_id, edges, nodes, data_edges_dir, label_exits, nodes_demand)
        feas_sols_clean[id_sol, :] = feas_sols_clean_id[:]
    return feas_sols_clean
