import numpy as np
from main_k_shortest_paths import *

def get_feas_sols_yens_exit_i(nodes, nodes_exit, nodes_exit_i, data_edges_dir, tij_edges_dir, demand_node, n_samples):

    data_edges_undir = np.vstack([data_edges_dir, data_edges_dir])
    data_edges_undir[(data_edges_dir.shape[0]):, 0] = data_edges_dir[:,1]
    data_edges_undir[(data_edges_dir.shape[0]):, 1] = data_edges_dir[:,0]
    n_k = n_samples
    tij_edges_dir = np.float32(tij_edges_dir)
    e_mat, e_mat_adj, t_mat, t_mat_adj = pre_graph(nodes, data_edges_dir, tij_edges_dir)
    e_mat_copy_all, e_mat_adj_copy_all, t_mat_copy_all, t_mat_adj_copy_all = pre_graph(nodes, data_edges_dir, tij_edges_dir)
    a_path_fin = []
    a_dis_exit_fin = []

    tij_edges_dir_copy = tij_edges_dir.copy()

    for e_i, exit_node in enumerate(nodes_exit_i):

        exit_node_rem = nodes_exit[nodes_exit != exit_node]
        
        chk1_i = np.isin(data_edges_dir[:,0], exit_node_rem)
        chk2_i = np.isin(data_edges_dir[:,1], exit_node_rem)
        chk12_i = (chk1_i | chk2_i)
        data_edges_dir_12 = data_edges_dir[~chk12_i]
        tij_edges_dir_12 = tij_edges_dir[~chk12_i]
        tij_edges_dir_copy_12 = tij_edges_dir_copy[~chk12_i]

        e_mat, e_mat_adj, t_mat, t_mat_adj = pre_graph(nodes, data_edges_dir_12, tij_edges_dir_copy_12)
        e_mat_copy_all, e_mat_adj_copy_all, t_mat_copy_all, t_mat_adj_copy_all = pre_graph(nodes, data_edges_dir_12, tij_edges_dir_copy_12)
        dis_fin, par, dis, cur_node, path_rev0, path0, dis_path0 = do_shortest(nodes, e_mat_adj, t_mat_adj, demand_node, exit_node)
        print(e_i, nodes_exit.shape)
        a_dis = [dis_path0]
        a_dis_exit = [dis_fin[cur_node]]
        a_path = [path0]

        b_dis = []
        b_dis_exit = []
        b_path = []

        do_l = True

        k_i = 0
        while do_l:
              k_i += 1
              #for k_i in range(1, n_k):
    
              a_dis_pre = a_dis[-1]
              a_path_pre = a_path[-1]
              for k_j in range(0, len(a_path_pre) - 1):
                 
                 demand_node_ki_j = a_path_pre[k_j]
                 dis_pre_ki_j = a_dis_pre[k_j]
                 path_pre_head = a_path_pre[:k_j]
                 dis_pre_head = a_dis_pre[:k_j]

                 path_chk_dem_nodes = []
                 path_chk = path_pre_head + [demand_node_ki_j]         
                 for k_i_pre in range(0, k_i):
                     path_chk_pre = a_path[k_i_pre][:len(path_chk)]
                     if (path_chk_pre == path_chk):
                        path_chk_dem_nodes.append(a_path[k_i_pre][len(path_chk)])
                 path_chk_dem_nodes = path_chk_dem_nodes + [a_path_pre[k_j+1]]

                 l1_chk = []
                 l1_node_chk = []

                 for l1, l1_node in enumerate(e_mat_adj[demand_node_ki_j]):
                     for l1_j, l1_node_j in enumerate(path_chk_dem_nodes):
                         if (l1_node == l1_node_j):
                            l1_chk.append(l1)
                            l1_node_chk.append(l1_node)

                 e_mat, e_mat_adj, t_mat, t_mat_adj = pre_graph(nodes, data_edges_dir_12, tij_edges_dir_copy_12)
                 demand_node_ki_j = a_path_pre[k_j]
                 dis_pre_ki_j = a_dis_pre[k_j]
                 path_pre_head = a_path_pre[:k_j]
                 dis_pre_head = a_dis_pre[:k_j]

                 for si1, snode1 in enumerate(path_pre_head):
                     for si2, snode2 in enumerate(e_mat_adj[snode1]):
                         t_mat_adj[snode1][si2] = np.inf

                 for l1_i in range(len(l1_chk)):
                     t_mat_adj[demand_node_ki_j][l1_chk[l1_i]] = np.inf

                 dis_fin_k_ij, par_k_ij, dis_k_ij, cur_node_k_ij, path_rev0_k_ij, path0_k_ij, dis_path_k_ij = do_shortest(nodes, e_mat_adj, t_mat_adj, demand_node_ki_j, exit_node)
                 if (dis_fin_k_ij[exit_node] != np.inf):
                     dis_prop = np.hstack([dis_pre_head, dis_path_k_ij + dis_pre_ki_j])
                     dis_exit_prop = dis_fin_k_ij[exit_node] + dis_pre_ki_j
                     path_prop = path_pre_head + path0_k_ij
           
                     chk1 = True
                     for pr_i in range(len(b_path)):
                         if (b_path[pr_i] == path_prop):
                            chk1 = False
                     if (chk1):
                          b_dis.append(dis_prop)
                          b_dis_exit.append(dis_exit_prop)
                          b_path.append(path_prop)

                 for si1, snode1 in enumerate(path_pre_head):
                     for si2, snode2 in enumerate(e_mat_adj[snode1]):
                         t_mat_adj[snode1][si2] = t_mat_adj_copy_all[snode1][si2]

              if (len(np.array(b_dis_exit)) > 0):
                 ag1 = np.argmin(np.array(b_dis_exit))
                 if (b_dis_exit[ag1] != np.inf):
                    a_path.append(b_path[ag1])
                    a_dis.append(b_dis[ag1])
                    a_dis_exit.append(b_dis_exit[ag1])
                 if (b_dis_exit[ag1] == np.inf):
                    do_l = False 
                 b_dis_exit[ag1] = np.inf
              else:
                 do_l = False
              if (k_i == n_k - 1):
                 do_l = False
              print(k_i, n_k-1, do_l, e_i, nodes_exit.shape)
        for fi1 in range(len(a_path)):
            a_path_fin.append(a_path[fi1])
            a_dis_exit_fin.append(a_dis_exit[fi1])
    a_dis_exit_fin_arr = np.array(a_dis_exit_fin)
    a_dis_exit_fin_argsort = np.argsort(a_dis_exit_fin_arr) #[0:n_k] (Note: uncomment #[0:n_k] to generate only n_k shortest)
    a_path_sorted = []
    for fi1 in range(0, len(a_dis_exit_fin_argsort)):
        a_path_sorted.append(a_path_fin[a_dis_exit_fin_argsort[fi1]])
    feas_sols_sorted = np.zeros((len(a_path_sorted), data_edges_undir.shape[0]), dtype=np.int32)
    for fi1 in range(0, len(a_dis_exit_fin_argsort)):
        fpath_i = a_path_sorted[fi1]
        for fi2 in range(0, len(fpath_i) - 1):
            for fi3 in range(data_edges_undir.shape[0]):
                if ((data_edges_undir[fi3][0] == fpath_i[fi2]) & (data_edges_undir[fi3][1] == fpath_i[fi2 + 1])):
                   feas_sols_sorted[fi1][fi3] = 1
    return feas_sols_sorted, a_path_sorted, a_dis_exit_fin_arr[a_dis_exit_fin_argsort]

def get_feas_sols_yens_distinct_exit_i(nodes, nodes_exit, nodes_exit_i, data_edges_dir, tij_edges_dir, demand_node, n_samples):

    data_edges_undir = np.vstack([data_edges_dir, data_edges_dir])
    data_edges_undir[(data_edges_dir.shape[0]):, 0] = data_edges_dir[:,1]
    data_edges_undir[(data_edges_dir.shape[0]):, 1] = data_edges_dir[:,0]
    n_k = n_samples
    tij_edges_dir = np.float32(tij_edges_dir)
    e_mat, e_mat_adj, t_mat, t_mat_adj = pre_graph(nodes, data_edges_dir, tij_edges_dir)
    e_mat_copy_all, e_mat_adj_copy_all, t_mat_copy_all, t_mat_adj_copy_all = pre_graph(nodes, data_edges_dir, tij_edges_dir)
    a_path_fin = []
    a_dis_exit_fin = []

    tij_edges_dir_copy = tij_edges_dir.copy()   

    for e_i, exit_node in enumerate(nodes_exit_i):

        exit_node_rem = nodes_exit[nodes_exit != exit_node]
        chk1_i = np.isin(data_edges_dir[:,0], exit_node_rem)
        chk2_i = np.isin(data_edges_dir[:,1], exit_node_rem)
        chk12_i = (chk1_i | chk2_i)
        data_edges_dir_12 = data_edges_dir[~chk12_i]
        tij_edges_dir_12 = tij_edges_dir_copy[~chk12_i]
        e_mat, e_mat_adj, t_mat, t_mat_adj = pre_graph(nodes, data_edges_dir_12, tij_edges_dir_12)
        e_mat_copy_all, e_mat_adj_copy_all, t_mat_copy_all, t_mat_adj_copy_all = pre_graph(nodes, data_edges_dir_12, tij_edges_dir_12)
        dis_fin, par, dis, cur_node, path_rev0, path0, dis_path0 = do_shortest(nodes, e_mat_adj, t_mat_adj, demand_node, exit_node)
        a_dis = [dis_path0]
        a_dis_exit = [dis_fin[cur_node]]
        a_path = [path0]

        b_dis = []
        b_dis_exit = []
        b_path = []

        do_l = True

        k_i = 0
        while do_l:
              k_i += 1
              path_k_i = a_path[-1]
              for fi2_ki in range(0, len(path_k_i) - 1):
                  for fi3_ki in range(data_edges_undir.shape[0]):
                      if ((data_edges_undir[fi3_ki][0] == path_k_i[fi2_ki]) & (data_edges_undir[fi3_ki][1] == path_k_i[fi2_ki + 1])):
                         t_mat_copy_all[data_edges_undir[fi3_ki][0], data_edges_undir[fi3_ki][1]] = (t_mat_copy_all[data_edges_undir[fi3_ki][0], data_edges_undir[fi3_ki][1]])*2.0
              e_mat_adj_copy_all, t_mat_adj_copy_all = pre_graph_only_tmat(nodes, e_mat_copy_all, t_mat_copy_all)
              dis_fin, par, dis, cur_node, path_rev0, path0, dis_path0 = do_shortest(nodes, e_mat_adj, t_mat_adj_copy_all, demand_node, exit_node)
              a_path.append(path0)
              a_dis_exit.append(dis_fin[cur_node])
              if (k_i == n_k - 1):
                 do_l = False
        for fi1 in range(len(a_path)):
            a_path_fin.append(a_path[fi1])
            a_dis_exit_fin.append(a_dis_exit[fi1])
    a_path_sorted = []
    a_dis_sorted = []
    for fi1 in range(0, len(a_path_fin)):
        a_path_sorted.append(a_path_fin[fi1])
        a_dis_sorted.append(a_dis_exit_fin[fi1])
    feas_sols_sorted = np.zeros((len(a_path_sorted), data_edges_undir.shape[0]), dtype=np.int32)
    for fi1 in range(0, len(a_path_sorted)):
        fpath_i = a_path_sorted[fi1]
        for fi2 in range(0, len(fpath_i) - 1):
            for fi3 in range(data_edges_undir.shape[0]):
                if ((data_edges_undir[fi3][0] == fpath_i[fi2]) & (data_edges_undir[fi3][1] == fpath_i[fi2 + 1])):
                   feas_sols_sorted[fi1][fi3] = 1
    return feas_sols_sorted, a_path_sorted, a_dis_sorted
