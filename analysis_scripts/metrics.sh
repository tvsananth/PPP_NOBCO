#!/bin/bash

RUN_DIR=$(pwd)
cd ${RUN_DIR}

ROOT_DIR=$(realpath "${RUN_DIR}/..")

module load anaconda3

python3 ./benchmark_statistics.py --path_results_dir_cinput=${ROOT_DIR}'/results' --json_filename_cinput='benchmark_10_0.75.json' --graph_size_cinput=10 --graph_den_cinput=0.75 --metric='TTS' --T_bnb_cinput=4320

python3 ./benchmark_statistics.py --path_results_dir_cinput=${ROOT_DIR}'/results' --json_filename_cinput='benchmark_10_0.75.json' --graph_size_cinput=10 --graph_den_cinput=0.75 --metric='ETS' --T_bnb_cinput=4320

python3 ./benchmark_statistics.py --path_results_dir_cinput=${ROOT_DIR}'/results' --json_filename_cinput='benchmark_20_0.75_hierch_decomp.json' --graph_size_cinput=20 --graph_den_cinput=0.75 --decomp_str_cinput='hierch_decomp' --metric='TTS' --T_bnb_cinput=72000

python3 ./benchmark_statistics.py --path_results_dir_cinput=${ROOT_DIR}'/results' --json_filename_cinput='benchmark_20_0.75_hierch_decomp.json' --graph_size_cinput=20 --graph_den_cinput=0.75 --decomp_str_cinput='hierch_decomp' --metric='ETS' --T_bnb_cinput=72000

python3 ./benchmark_statistics.py --path_results_dir_cinput=${ROOT_DIR}'/results' --json_filename_cinput='benchmark_20_0.75_alt_min.json' --graph_size_cinput=20 --graph_den_cinput=0.75 --decomp_str_cinput='alt_min' --metric='TTS' --T_bnb_cinput=72000

python3 ./benchmark_statistics.py --path_results_dir_cinput=${ROOT_DIR}'/results' --json_filename_cinput='benchmark_20_0.75_alt_min.json' --graph_size_cinput=20 --graph_den_cinput=0.75 --decomp_str_cinput='alt_min' --metric='ETS' --T_bnb_cinput=72000
