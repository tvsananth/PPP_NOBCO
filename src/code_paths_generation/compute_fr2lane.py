import numpy as np
from sim_anneal_paths import *
from compute_A_undirected import *
from compute_sols_undir import *
from clean_cycles_sols import *
from graver_diff_sols import *
from edges_cost import *
import time
import argparse

def load_data_graph(save_dir):

    data_nodes = np.genfromtxt(save_dir + '/data_nodes.txt')
    data_edges_dir = np.int32(np.genfromtxt(save_dir + '/data_edges_dir.txt'))
    cij_edges_dir = np.genfromtxt(save_dir + '/cij_edges_dir.txt')
    tij_edges_dir = np.genfromtxt(save_dir + '/tij_edges_dir.txt')
    lanes_edges_dir = np.int32(np.genfromtxt(save_dir + '/lanes_edges_dir.txt'))
    nodes_demand = np.int32(np.genfromtxt(save_dir + '/nodes_demand.txt'))
    nodes_exit = np.int32(np.genfromtxt(save_dir + '/nodes_exit.txt'))

    return data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit


def do_undir(cij_edges_dir, tij_edges_dir, lanes_edges_dir):

    cij_edges_undir = np.hstack([cij_edges_dir, cij_edges_dir])

    tij_edges_undir = np.hstack([tij_edges_dir, tij_edges_dir])

    lanes_edges_undir = np.hstack([lanes_edges_dir, lanes_edges_dir])

    return cij_edges_undir, tij_edges_undir, lanes_edges_undir

