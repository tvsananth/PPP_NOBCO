import numpy as np
import csv
import pandas as pd

def get_files(input_dir, graph_size='10', graph_den='0.25', graph_instance_index='0'):
    nodes_file = input_dir + '/random_instances/nodes_' + graph_size + '_' + graph_den + '_' + graph_instance_index + '.csv'
    edges_file = input_dir + '/random_instances/edges_' + graph_size + '_' + graph_den + '_' + graph_instance_index + '.csv'

    return nodes_file, edges_file


def read_nodes_random(nodes_file):

    read_nodes = []
    with open(nodes_file, 'r') as nodes_csv_file:
         reader = csv.reader(nodes_csv_file)
         for r in reader:
             read_nodes.append(r)
             print(r)
         nodes_csv_file.close()
    n_nodes = len(read_nodes[1:])

    return read_nodes, n_nodes


def read_edges_random(edges_file):
    data_edges = []
    with open(edges_file, 'r') as edges_csv_file:
         reader = csv.reader(edges_csv_file)
         for r in reader:
             data_edges.append(r)
             print(r)
         edges_csv_file.close()
    n_edges = len(data_edges[1:])//2
    data_edges_re = []
    data_edges_re.append(data_edges[0])
    for i1 in range(n_edges):
        data_edges_re.append(data_edges[1+i1*2])
    for i1 in range(n_edges):
        data_edges_re.append(data_edges[1+i1*2 + 1])

    return data_edges, data_edges_re, n_edges


def clean_data_random(graph_size, graph_den, graph_instance_index, input_dir):

    nodes_file, edges_file = get_files(input_dir, graph_size=graph_size, graph_den=graph_den, graph_instance_index=graph_instance_index)

    read_nodes, n_nodes = read_nodes_random(nodes_file)
    data_edges, data_edges_re, n_edges = read_edges_random(edges_file)
    
    df_nodes = pd.DataFrame(read_nodes[1:])
    df_nodes.to_csv(f'{input_dir}/random_instances_reorder/nodes_{graph_size}_{graph_den}_{graph_instance_index}.csv', header=['x','y','SR_demand','exits','injured'], index=False)

    df_edges = pd.DataFrame(data_edges_re[1:])
    df_edges.to_csv(f'{input_dir}/random_instances_reorder/edges_{graph_size}_{graph_den}_{graph_instance_index}.csv', header=['from','to','t0','c0','lanes','initial_SR','initial_FR'], index=False)

input_dir = '../data/instances'
"""
for i1, graph_size in enumerate(['10']):
    for j1, graph_den in enumerate(['0.75']):
        for k1, graph_instance_index in enumerate(np.arange(10)):
            clean_data_random(graph_size, graph_den, str(graph_instance_index), input_dir)
"""
