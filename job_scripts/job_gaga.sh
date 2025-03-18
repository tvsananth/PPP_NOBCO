#!/bin/bash

RUN_DIR=$(pwd)
cd ${RUN_DIR}

INSTANCE_NSIZE=10
INSTANCE_ESIZE=0.75
INSTANCE_NUM=0
for TASK_ID in $(seq 0 99)
do
	bash run_gaga.sh ${INSTANCE_NSIZE} ${INSTANCE_ESIZE} ${INSTANCE_NUM} ${TASK_ID} &
        sleep 1	
done

wait
