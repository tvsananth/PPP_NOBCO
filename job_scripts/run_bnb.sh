#!/bin/bash

INSTANCE_NSIZE=$1
INSTANCE_ESIZE=$2
RUN_TIME=$3
INSTANCE_NUM=$4

ISTR=${INSTANCE_NSIZE}_${INSTANCE_ESIZE}_${INSTANCE_NUM}
GSTR=graph_${ISTR}

RUN_DIR=$(pwd)
ROOT_DIR=$(realpath "${RUN_DIR}/..")
SRC_DIR=${ROOT_DIR}/src
DATA_DIR=${ROOT_DIR}/data
LOGS_DIR=${ROOT_DIR}/logs
SCRIPTS_DIR=${ROOT_DIR}/scripts
INSTANCES_DIR=${DATA_DIR}/instances/random_instances
INSTANCES_REORDER_DIR=${DATA_DIR}/instances/random_instances_reorder
OUTPUT_DIR=${DATA_DIR}/output/random_instances

cd ${SRC_DIR}/code_bnb

pwd

timeout ${RUN_TIME} ./leblanc_solver ${INSTANCES_DIR}/nodes_${ISTR}.csv ${INSTANCES_DIR}/edges_${ISTR}.csv FR ${OUTPUT_DIR}/branch_and_bound/logs_${ISTR}.npy > ${LOGS_DIR}/logs_bnb/logs_${ISTR}.csv 
