#!/bin/bash

#PROJECT_HOME="/Users/joefresna/DecisionOnNetworks"
PROJECT_HOME="/home/ac1ar/DecisionsOnNetworks"
CONF_DIR="${PROJECT_HOME}/conf_cluster"
mkdir -p ${CONF_DIR}
TEMPLATE_SETTINGS="${PROJECT_HOME}/conf/DecNet.template.config"
PYPATH="${PROJECT_HOME}/src/"
EXEC_FILE="${PROJECT_HOME}/src/DecNet/DecisionProcess.py"
OUTPUT_DATA_DIR="${PROJECT_HOME}/conf1/"

SEED=2289
NUM_EXP=100
NUM_RUN=100

AGENT_TYPE='simple'
DEC_MODEL_LIST=('conf-perfect')

### Unnecessary for this analysis 
DDM_DRIFT_DIST='from_accuracy'
DDM_BASE_DRIFT=0.1
DDM_DRIFT_RANGE_MIN=-1
DDM_DRIFT_RANGE_MAX=1
DDM_DRIFT_STD_DEV=0.3
#THRESHOLD=1
DDM_NOISE_STD_DEV=0.5
INTERR_TIME=2

#NUM_NODES_LIST=(100 500 1000)
NUM_NODES_LIST=($(seq 11 4 31))
NET_TYPE='full'
LINK_PROBABILITY=1
NUM_EDGES=1

ACC_MEAN_LIST=( 0.6 )
ACC_STD_DEV_LIST=( 0.3 )

COUNT=0

for DEC_MODEL in ${DEC_MODEL_LIST[*]}
do
	for NUM_NODES in ${NUM_NODES_LIST[*]}
	do
		for ACC_MEAN in ${ACC_MEAN_LIST[*]}
		do
			for ACC_STD_DEV in ${ACC_STD_DEV_LIST[*]}
			do
				
				JOB_PARAM="net-${NET_TYPE}_${AGENT_TYPE}_nodes-${NUM_NODES}_model-${DEC_MODEL}_acc-${ACC_MEAN}_acstdv-${ACC_STD_DEV}"
				OUT_TXT="${OUTPUT_DATA_DIR}out_${JOB_PARAM}.txt"
				#OUT_PDF="${OUTPUT_DATA_DIR}pdf/out_${JOB_PARAM}"
					
				CONF_FILE="${CONF_DIR}/decnet_${JOB_PARAM}.config"
					
				sed -e "s|SEED|${SEED}|" \
					-e "s|NUM_EXP|${NUM_EXP}|" \
					-e "s|NUM_RUN|${NUM_RUN}|" \
					-e "s|AGENT_TYPE|${AGENT_TYPE}|" \
					-e "s|ACC_MEAN|${ACC_MEAN}|" \
					-e "s|ACC_STD_DEV|${ACC_STD_DEV}|" \
					-e "s|DEC_MODEL|${DEC_MODEL}|" \
					-e "s|DDM_DRIFT_DIST|${DDM_DRIFT_DIST}|" \
					-e "s|DDM_BASE_DRIFT|${DDM_BASE_DRIFT}|" \
					-e "s|DDM_DRIFT_RANGE_MIN|${DDM_DRIFT_RANGE_MIN}|" \
					-e "s|DDM_DRIFT_RANGE_MAX|${DDM_DRIFT_RANGE_MAX}|" \
					-e "s|DDM_DRIFT_STD_DEV|${DDM_DRIFT_STD_DEV}|" \
					-e "s|DDM_NOISE_STD_DEV|${DDM_NOISE_STD_DEV}|" \
					-e "s|THRESHOLD|${THRESHOLD}|" \
					-e "s|INTERR_TIME|${INTERR_TIME}|" \
					-e "s|NET_TYPE|${NET_TYPE}|" \
					-e "s|NUM_NODES|${NUM_NODES}|" \
					-e "s|NUM_EDGES|${NUM_EDGES}|" \
					-e "s|LINK_PROBABILITY|${LINK_PROBABILITY}|" \
					-e "s|OUT_TXT|${OUT_TXT}|" \
						${TEMPLATE_SETTINGS} > ${CONF_FILE}
							
				export PYTHONPATH=${PYPATH}
				#COMMAND="python3 ${EXEC_FILE} ${CONF_FILE}"
				#COMMAND="./run_job.sh ${EXEC_FILE} ${CONF_FILE}"
				COMMAND="qsub run_job.sh ${EXEC_FILE} ${CONF_FILE}"
				${COMMAND}
				COUNT=$((COUNT + 1))
			done
		done
	done
done

echo "Submitted " $COUNT " jobs"
