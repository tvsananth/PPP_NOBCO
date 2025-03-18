import numpy as np
import argparse
import json

def main_compute():

    parser = argparse.ArgumentParser()

    parser.add_argument('-path_output_data_dir', '--path_output_data_dir_cinput', type=str, help='path to the Output data directory')
    parser.add_argument('-graph_size', '--graph_size_cinput', type=int, default=0, help='graph node size')
    parser.add_argument('-graph_den', '--graph_den_cinput', type=float, default=0, help='graph edge density')
    parser.add_argument('-decomp_str', '--decomp_str_cinput', type=str, default='', help='Type of decomposition used')
    parser.add_argument('-path_results_dir', '--path_results_dir_cinput', type=str, default='', help='directory to save JSON files')

    args = parser.parse_args()
    path_output_data_dir = args.path_output_data_dir_cinput
    graph_size = args.graph_size_cinput
    graph_den  = args.graph_den_cinput
    decomp_str = args.decomp_str_cinput
    path_results_dir = args.path_results_dir_cinput

    ini = 0
    fin = 99
    
    path_save_data_prec  = path_output_data_dir + '/random_instances/graph_' + str(graph_size) + '_' + str(graph_den) + '_'

    if (graph_size == 20):
        decomp_str = '_' + decomp_str

    solv_params = {}
    solv_params['M_seeds'] = 100

    ins_dict = []
    ins_list = np.arange(0,10)
    for i_ins_id_num, ins_id in enumerate(ins_list):

        bnb_a1 = np.fromfile(path_output_data_dir + '/random_instances/branch_and_bound/logs_' + str(graph_size) + '_' + str(graph_den) + '_' + str(ins_id) + '.npy')
        bnb_size = bnb_a1.shape[0]//3
        bnb_data = bnb_a1[0:bnb_size*3].reshape((bnb_size, 3))
        bnb_time = bnb_data[-1][2]
        bnb_obj  = bnb_data[-1][1]
        path_save_data = path_save_data_prec + str(ins_id) + '_run'
        print(path_save_data, ini, fin, bnb_obj, bnb_time)
   
        tts_sols = np.zeros(fin - ini + 1)
        ets_sols = np.zeros(fin - ini + 1)
        fr_sols  = np.zeros(fin - ini + 1)
        for i_id in range(ini, fin+1):
            out_dir = path_save_data + str(i_id) + decomp_str
            t_sa_g = np.genfromtxt(out_dir + '/times_sa_graver.txt')
            t_sa = t_sa_g[0]
            t_g  = t_sa_g[1]
            if (graph_size == 20):
                t_analog = np.genfromtxt(out_dir + '/times_paths_only_analog.txt')
            else:
                t_analog = t_sa
            t_gaga = np.fromfile(out_dir + '/Graver_walk/e21time_sol_srallFR_op1000.npy')
            t_leb  = np.fromfile(out_dir + '/Graver_walk/e21sr_time_nogama_op1000.npy')
            fr_sol = np.fromfile(out_dir + '/Graver_walk/e21sr_sol_nogama_op1000.npy')
            chk_sol = (fr_sol > 0)
            t_gaga = t_gaga[chk_sol]
            t_leb  = t_leb[chk_sol]
            fr_sol = fr_sol[chk_sol]
            tts_id = t_sa + t_g + (t_gaga.sum()) + (t_leb.sum())
            ets_id = (t_sa - t_analog) + t_g + (t_gaga.sum()) + (t_leb.sum())
            sol_id = fr_sol.min()
            tts_sols[i_id] = tts_id
            ets_sols[i_id] = ets_id
            fr_sols[i_id]  = sol_id
            #print(i_id, tts_id, ets_id, sol_id)
        fr_sols_diff = (bnb_obj - fr_sols)*100.0/bnb_obj 
        fr_sols_solved = (fr_sols <= bnb_obj).astype(int)
        n_sols_solved = fr_sols_solved.sum()
        ins_dict_id = {}
        ins_dict_id['solver'] = 'GAGA'
        ins_dict_id['solver parameters'] = solv_params
        ins_dict_id['hardware'] = 'cpu:intel'
        ins_dict_id['set'] = 'RandomGraphs: n = ' + str(graph_size) + ', p = ' + str(graph_den) 
        ins_dict_id['instance_idx'] = str(ins_id)
        ins_dict_id['cutoff_type'] = 'Algorithm Completed'
        ins_dict_id['cutoff'] = 'Completed'
        ins_dict_id['runs_attempted'] = 100
        ins_dict_id['runs_solved'] = str(n_sols_solved)
        ins_dict_id['n_unsat_clauses'] = (1 - fr_sols_solved).tolist()
        ins_dict_id['Objective function'] = (fr_sols).tolist()
        ins_dict_id['cSOA Objective function'] = (bnb_obj).tolist()
        ins_dict_id['configurations'] = fr_sols_diff.tolist()
        ins_dict_id['pre_runtime_seconds'] = tts_sols.tolist()
        ins_dict_id['runtime_seconds'] = tts_sols.tolist()
        ins_dict_id['energy'] = ets_sols.tolist()
        ins_dict.append(ins_dict_id) 
        print(ins_id)         
    with open(path_results_dir + '/benchmark_' + str(graph_size) + '_' + str(graph_den) + decomp_str + '.json', 'w') as f_out:
        json.dump(ins_dict, f_out)
    with open(path_results_dir + '/benchmark_' + str(graph_size) + '_' + str(graph_den) + decomp_str + '.json', 'r') as f_in:
        ins_data_jfile = json.load(f_in)
    
if __name__ == '__main__':
     main_compute()
