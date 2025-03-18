import numpy as np
import argparse
import json

def main_compute():

    parser = argparse.ArgumentParser()

    parser.add_argument('-path_output_data_dir', '--path_output_data_dir_cinput', type=str, help='path to the Output data directory')
    parser.add_argument('-graph_size', '--graph_size_cinput', type=int, default=0, help='graph node size')
    parser.add_argument('-graph_den', '--graph_den_cinput', type=float, default=0, help='graph edge density')
    parser.add_argument('-path_results_dir', '--path_results_dir_cinput', type=str, default='', help='directory to save JSON files')
    parser.add_argument('-run_time', '--run_time_cinput', type=float, default=0, help='Run time of branch-and-bound')

    args = parser.parse_args()
    path_output_data_dir = args.path_output_data_dir_cinput
    graph_size = args.graph_size_cinput
    graph_den  = args.graph_den_cinput
    path_results_dir = args.path_results_dir_cinput
    run_time = args.run_time_cinput

    solv_params = {}
    ins_dict = []
    ins_list = np.arange(0,10)
    for i_ins_id_num, ins_id in enumerate(ins_list):
        bnb_a1 = np.fromfile(path_output_data_dir + '/random_instances/branch_and_bound/logs_' + str(graph_size) + '_' + str(graph_den) + '_' + str(ins_id) + '.npy')
        bnb_size = bnb_a1.shape[0]//3
        bnb_data = bnb_a1[0:bnb_size*3].reshape((bnb_size, 3))
        bnb_time = bnb_data[-1][2]
        bnb_obj  = bnb_data[-1][1]
        print(i_ins_id_num, bnb_obj, bnb_time)

        ins_dict_id = {}
        ins_dict_id['solver'] = 'Branch-and-Bound'
        ins_dict_id['solver parameters'] = solv_params
        ins_dict_id['hardware'] = 'cpu:intel'
        ins_dict_id['set'] = 'RandomGraphs: n = ' + str(graph_size) + ', p = ' + str(graph_den)
        ins_dict_id['instance_idx'] = str(ins_id)
        ins_dict_id['cutoff_type'] = 'Run time'
        ins_dict_id['cutoff'] = str(run_time) + ' seconds'
        ins_dict_id['runs_attempted'] = 1
        ins_dict_id['cSOA objectve function'] = bnb_obj
        ins_dict_id['runtime_seconds'] = run_time
        ins_dict_id['energy'] = run_time
        ins_dict.append(ins_dict_id) 
        print(ins_id)         
    with open(path_results_dir + '/benchmark_' + str(graph_size) + '_' + str(graph_den) + '_bnb.json', 'w') as f_out:
        json.dump(ins_dict, f_out)
    with open(path_results_dir + '/benchmark_' + str(graph_size) + '_' + str(graph_den) + '_bnb.json', 'r') as f_in:
        ins_data_jfile = json.load(f_in)
    
if __name__ == '__main__':
     main_compute()
