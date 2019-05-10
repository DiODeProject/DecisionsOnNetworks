#!/bin/bash

#PROJECT_HOME="/Users/joefresna/DecisionOnNetworks"
PROJECT_HOME="/home/ac1ar/DecisionsOnNetworks"
CONF_DIR="${PROJECT_HOME}/conf_cluster"
mkdir -p ${CONF_DIR}
TEMPLATE_SETTINGS="${PROJECT_HOME}/conf/DecNet.template.config"
PYPATH="${PROJECT_HOME}/src/"
EXEC_FILE="${PROJECT_HOME}/src/DecNet/DecisionProcess.py"
OUTPUT_DATA_DIR="${PROJECT_HOME}/belief8/"

SEED=2289
NUM_EXP=10000
NUM_RUN=1
MAX_TIME=200

AGENT_TYPE='simple'
DEC_MODEL_LIST=('belief-init' 'conf-perfect')
UPDATE_CONF_LIST=('belief-up' 'optim-up')
#BELIEF_EPSILON_LIST=($(seq 0.07 0.001 0.1))
BELIEF_EPSILON_LIST=($(seq 0.01 0.005 0.15))

### Unnecessary for this analysis 
DDM_DRIFT_DIST='normal'
DDM_BASE_DRIFT=0.1
DDM_DRIFT_RANGE_MIN=-1
DDM_DRIFT_RANGE_MAX=1
DDM_DRIFT_STD_DEV=0.3
DDM_NOISE_STD_DEV=0.5
#THRESHOLD=1
INTERR_TIME=10

#NUM_NODES=20
#NUM_NODES_LIST=(31)
NUM_NODES_LIST=($(seq 11 12 47))
NET_TYPE_LIST=('full' 'erdos-renyi' 'barabasi-albert')
#NET_TYPE_LIST=('barabasi-albert')
#NET_TYPE_LIST=('erdos-renyi')
#LINK_PROBABILITY_LIST=($(seq 0.3 0.1 0.8))
#NUM_EDGES_LIST=($(seq 2 2 10))
LINK_PROBABILITY_LIST=(0.2 0.4 0.6 0.8)
NUM_EDGES_LIST=(3 7 11 15)
#LINK_PROBABILITY_LIST=(0.3)
#NUM_EDGES_LIST=(16)

ACC_MEAN_LIST=( 0.6 )
ACC_STD_DEV_LIST=( 0.12 )
ACC_TRUNCATED='true'

COUNT=0

for DEC_MODEL in ${DEC_MODEL_LIST[*]}
do
	for NET_TYPE in ${NET_TYPE_LIST[*]}
	do
		for NUM_NODES in ${NUM_NODES_LIST[*]}
		do
			for NUM_EDGES in ${NUM_EDGES_LIST[*]}
			do
				for LINK_PROBABILITY in ${LINK_PROBABILITY_LIST[*]}
				do
					for ACC_MEAN in ${ACC_MEAN_LIST[*]}
					do
						for ACC_STD_DEV in ${ACC_STD_DEV_LIST[*]}
						do
							for UPDATE_CONF in ${UPDATE_CONF_LIST[*]}
							do
								for BELIEF_EPSILON in ${BELIEF_EPSILON_LIST[*]}
								do
									### Skipping unnecessary iterations
									if [ "$NET_TYPE" != "erdos-renyi" ] && [ $LINK_PROBABILITY != ${LINK_PROBABILITY_LIST[0]} ];
									then
										continue
									fi
									if [ "$NET_TYPE" != "barabasi-albert" ] && [ $NUM_EDGES != ${NUM_EDGES_LIST[0]} ];
									then
										continue
									fi
									if [ "$DEC_MODEL" == "belief-init" ] && [ $UPDATE_CONF != "belief-up" ];
									then
										continue
									fi
									if [ "$UPDATE_CONF" != "belief-up" ] && [ $BELIEF_EPSILON != ${BELIEF_EPSILON_LIST[0]} ];
									then
										continue
									fi
									
									
									JOB_PARAM="net-${NET_TYPE}_nodes-${NUM_NODES}_edges-${NUM_EDGES}_link-${LINK_PROBABILITY}_model-${DEC_MODEL}_up-${UPDATE_CONF}_acc-${ACC_MEAN}_acstdv-${ACC_STD_DEV}_eps-${BELIEF_EPSILON}"
									OUT_TXT="${OUTPUT_DATA_DIR}out_${JOB_PARAM}.txt"
										
									CONF_FILE="${CONF_DIR}/decnet_${JOB_PARAM}.config"
										
									sed -e "s|SEED|${SEED}|" \
										-e "s|NUM_EXP|${NUM_EXP}|" \
										-e "s|NUM_RUN|${NUM_RUN}|" \
										-e "s|MAX_TIME|${MAX_TIME}|" \
										-e "s|AGENT_TYPE|${AGENT_TYPE}|" \
										-e "s|ACC_MEAN|${ACC_MEAN}|" \
										-e "s|ACC_STD_DEV|${ACC_STD_DEV}|" \
										-e "s|ACC_TRUNCATED|${ACC_TRUNCATED}|" \
										-e "s|DEC_MODEL|${DEC_MODEL}|" \
										-e "s|UPDATE_CONF|${UPDATE_CONF}|" \
										-e "s|BELIEF_EPSILON|${BELIEF_EPSILON}|" \
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
									#${COMMAND}
									while ! ${COMMAND}
									do
										sleep 2
									done
										
									COUNT=$((COUNT + 1))
								done
							done
						done
					done
				done
			done
		done
	done
done

echo "Submitted " $COUNT " jobs"
