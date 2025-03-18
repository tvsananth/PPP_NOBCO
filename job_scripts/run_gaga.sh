#!/bin/bash

INSTANCE_NSIZE=$1
INSTANCE_ESIZE=$2
INSTANCE_NUM=$3
TASK_ID=$4

ISTR=${INSTANCE_NSIZE}_${INSTANCE_ESIZE}_${INSTANCE_NUM}
GSTR=graph_${ISTR}

RUN_DIR=$(pwd)
ROOT_DIR=$(realpath "${RUN_DIR}/..")
SRC_DIR=${ROOT_DIR}/src
DATA_DIR=${ROOT_DIR}/data
LOGS_DIR=${ROOT_DIR}/logs
INSTANCES_DIR=${DATA_DIR}/instances/random_instances
INSTANCES_REORDER_DIR=${DATA_DIR}/instances/random_instances_reorder
OUTPUT_DIR=${DATA_DIR}/output/random_instances

cd ${SRC_DIR}/code_paths_generation

pwd

python3 ./main_gen_paths_random_instances.py --nodes_file_cinput=${INSTANCES_DIR}'/nodes_'${ISTR}'.csv' --edges_file_cinput=${INSTANCES_DIR}'/edges_'${ISTR}'.csv' --path_save_data_cinput=${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID} --n_samples_cinput=10000 --feas_sols_cinput=True --graver_basis_cinput=True > ${LOGS_DIR}/logs_gaga/out_paths_${ISTR}_${TASK_ID}.log

cd ${SRC_DIR}/code_random

pwd

./main_FR_random_instances ${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID} ${INSTANCES_REORDER_DIR}'/nodes_'${ISTR}'.csv' ${INSTANCES_REORDER_DIR}'/edges_'${ISTR}'.csv' > ${LOGS_DIR}/logs_gaga/out_augment_${ISTR}_${TASK_ID}.log

cd ${RUN_DIR}

ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/data_nodes.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/nodes_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/*edges*dir.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/feas_sols_sorted_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/fsol_list_reduced_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/fsol_graver_r.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/ini_fr_rand_testfr.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/ini_sol_rand_testfr.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/rand_ini_fr_dir
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/rand_ini_sols_dir
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/Graver_walk/FR_*.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/Graver_walk/SR_*.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/Graver_walk/cijSR_*.npy
	
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/data_nodes.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/nodes_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/*edges*dir.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/feas_sols_sorted_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/fsol_list_reduced_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/fsol_graver_r.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/ini_fr_rand_testfr.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/ini_sol_rand_testfr.npy
rm -r ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/rand_ini_fr_dir
rm -r ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/rand_ini_sols_dir
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/Graver_walk/FR_*.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/Graver_walk/SR_*.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}/Graver_walk/cijSR_*.npy
