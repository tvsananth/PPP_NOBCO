#!/bin/bash

RUN_DIR=$(pwd)
cd ${RUN_DIR}

INSTANCE_NSIZE=20
INSTANCE_ESIZE=0.75
INSTANCE_NUM=0
DECOMP_STR='hierch_decomp'
##DECOMP_STR='alt_min', 'hierch_decomp'

for TASK_ID in $(seq 0 99)
do
	bash run_dgaga.sh ${INSTANCE_NSIZE} ${INSTANCE_ESIZE} ${INSTANCE_NUM} ${DECOMP_STR} ${TASK_ID} &
        sleep 1	
done

wait
