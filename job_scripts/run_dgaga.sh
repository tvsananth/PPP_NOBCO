#!/bin/bash

INSTANCE_NSIZE=$1
INSTANCE_ESIZE=$2
INSTANCE_NUM=$3
DECOMP_STR=$4
TASK_ID=$5

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

python3 ./main_gen_paths_decompose_random_instances_part_ini.py --nodes_file_cinput=${INSTANCES_DIR}'/nodes_'${ISTR}'.csv' --edges_file_cinput=${INSTANCES_DIR}'/edges_'${ISTR}'.csv' --path_save_data_cinput=${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID} --decompose_type_cinput=${DECOMP_STR} --n_samples_cinput=10000 --feas_sols_cinput=True --graver_basis_cinput=True > ${LOGS_DIR}/logs_dgaga/out_paths_ini_${ISTR}_${DECOMP_STR}_${TASK_ID}.log

if [[ "${DECOMP_STR}" == "alt_min" ]]
then
	python3 ./decompose_random_graph_alt_min_outward.py --nodes_file_cinput=${INSTANCES_DIR}'/nodes_'${ISTR}'.csv' --edges_file_cinput=${INSTANCES_DIR}'/edges_'${ISTR}'.csv' --path_save_data_cinput=${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID} --decompose_type_cinput=${DECOMP_STR} --max_cluster_size_cinput=40 > ${LOGS_DIR}/logs_dgaga/out_paths_decompose_${ISTR}_${DECOMP_STR}_${TASK_ID}.log
fi

if [[ "${DECOMP_STR}" == "hierch_decomp" ]]
then
	python3 ./decompose_random_graph_hierch_decomp.py --nodes_file_cinput=${INSTANCES_DIR}'/nodes_'${ISTR}'.csv' --edges_file_cinput=${INSTANCES_DIR}'/edges_'${ISTR}'.csv' --path_save_data_cinput=${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID} --decompose_type_cinput=${DECOMP_STR} --max_cluster_size_cinput=40  > ${LOGS_DIR}/logs_dgaga/out_paths_decompose_${ISTR}_${DECOMP_STR}_${TASK_ID}.log
fi

python3 ./main_gen_paths_decompose_random_instances_part_fin.py --nodes_file_cinput=${INSTANCES_DIR}'/nodes_'${ISTR}'.csv' --edges_file_cinput=${INSTANCES_DIR}'/edges_'${ISTR}'.csv' --path_save_data_cinput=${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID} --decompose_type_cinput=${DECOMP_STR} --n_samples_cinput=10000 --feas_sols_cinput=True --graver_basis_cinput=True > ${LOGS_DIR}/logs_dgaga/out_paths_fin_${ISTR}_${DECOMP_STR}_${TASK_ID}.log

cd ${SRC_DIR}/code_random
pwd
#echo ${SRC_DIR}/code_random

./main_FR_random_instances ${OUTPUT_DIR}'/graph_'${ISTR}'_run'${TASK_ID}'_'${DECOMP_STR} ${INSTANCES_REORDER_DIR}'/nodes_'${ISTR}'.csv' ${INSTANCES_REORDER_DIR}'/edges_'${ISTR}'.csv' > ${LOGS_DIR}/logs_dgaga/out_augment_${ISTR}_${DECOMP_STR}_${TASK_ID}.log

cd ${RUN_DIR}

ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/data_nodes.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/nodes_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/*edges*dir.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/rand_sols_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/feas_sols_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/fsol_list_reduced_*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/time_onlySA_par_cluster*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/times_onlySA_full_demand*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/times_par_full_demand*.txt
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/fsol_graver_r.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/ini_fr_rand_testfr.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/ini_sol_rand_testfr.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/rand_ini_fr_dir
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/rand_ini_sols_dir
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/Graver_walk/FR_*.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/Graver_walk/SR_*.npy
ls ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/Graver_walk/cijSR_*.npy

rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/data_nodes.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/nodes_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/*edges*dir.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/rand_sols_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/feas_sols_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/fsol_list_reduced_*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/time_onlySA_par_cluster*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/times_onlySA_full_demand*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/times_par_full_demand*.txt
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/fsol_graver_r.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/ini_fr_rand_testfr.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/ini_sol_rand_testfr.npy
rm -r ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/rand_ini_fr_dir
rm -r ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/rand_ini_sols_dir
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/Graver_walk/FR_*.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/Graver_walk/SR_*.npy
rm ${OUTPUT_DIR}/graph_${ISTR}_run${TASK_ID}_${DECOMP_STR}/Graver_walk/cijSR_*.npy

