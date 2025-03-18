#!/bin/bash

RUN_DIR=$(pwd)
cd ${RUN_DIR}

ROOT_DIR=$(realpath "${RUN_DIR}/..")

python3 ./make_json_files_gaga.py --path_output_data_dir_cinput=${ROOT_DIR}'/data/output' --graph_size_cinput=10 --graph_den_cinput=0.75 --path_results_dir_cinput=${ROOT_DIR}'/results'

python3 ./make_json_files_gaga.py --path_output_data_dir_cinput=${ROOT_DIR}'/data/output' --graph_size_cinput=20 --graph_den_cinput=0.75 --decomp_str_cinput='alt_min' --path_results_dir_cinput=${ROOT_DIR}'/results'

python3 ./make_json_files_gaga.py --path_output_data_dir_cinput=${ROOT_DIR}'/data/output' --graph_size_cinput=20 --graph_den_cinput=0.75 --decomp_str_cinput='hierch_decomp' --path_results_dir_cinput=${ROOT_DIR}'/results'

python3 ./make_json_files_bnb.py --path_output_data_dir_cinput=${ROOT_DIR}'/data/output' --graph_size_cinput=10 --graph_den_cinput=0.75 --path_results_dir_cinput=${ROOT_DIR}'/results' --run_time_cinput=4320

python3 ./make_json_files_bnb.py --path_output_data_dir_cinput=${ROOT_DIR}'/data/output' --graph_size_cinput=20 --graph_den_cinput=0.75 --path_results_dir_cinput=${ROOT_DIR}'/results' --run_time_cinput=72000

