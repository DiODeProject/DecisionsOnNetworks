#!/bin/bash

#PROJECT_HOME="/Users/joefresna/DecisionOnNetworks"
PROJECT_HOME="/home/ac1ar/DecisionsOnNetworks"
CONF_DIR="${PROJECT_HOME}/conf_cluster"
mkdir -p ${CONF_DIR}
TEMPLATE_SETTINGS="${PROJECT_HOME}/conf/DecNet.template.config"
PYPATH="${PROJECT_HOME}/src/"
EXEC_FILE="${PROJECT_HOME}/src/DecNet/DecisionProcess.py"
OUTPUT_DATA_DIR="${PROJECT_HOME}/nets-1708/"

SEED=2289
NUM_EXP=10000
NUM_RUN=1
MAX_TIME=200

AGENT_TYPE='simple'
#DEC_MODEL_LIST=('conf-perfect' 'majority-rand' 'majority-bias' 'majority-inhibit' 'best-acc')
#UPDATE_CONF_LIST=('no-up' 'theta-up' 'theta-norm' 'optim-up')
DEC_MODEL_LIST=('conf-perfect' 'majority-rand' 'best-acc')
#DEC_MODEL_LIST=('conf-perfect')
UPDATE_CONF_LIST=('no-up' 'optim-up' 'belief-up')
#UPDATE_CONF_LIST=('belief-up')

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
#NUM_NODES_LIST=(100 500 1000)
NUM_NODES_LIST=($(seq 11 6 47))
NET_TYPE_LIST=('full' 'erdos-renyi' 'barabasi-albert')
#NET_TYPE_LIST=('full')
#NET_TYPE_LIST=('erdos-renyi' 'barabasi-albert')
LINK_PROBABILITY_LIST=($(seq 0.2 0.1 0.8))
#LINK_PROBABILITY_LIST=(0.2)
#NUM_EDGES_LIST=(3)
NUM_EDGES_LIST=($(seq 3 2 15))

ACC_MEAN_LIST=( 0.6 )
ACC_STD_DEV_LIST=( 0.12 )
#ACC_TRUNCATED_LIST=('true' 'false')
ACC_TRUNCATED='true'

COUNT=0

for DEC_MODEL in ${DEC_MODEL_LIST[*]}
do
	for NET_TYPE in ${NET_TYPE_LIST[*]}
	do
		for NUM_NODES in ${NUM_NODES_LIST[*]}
		do
			BELIEF_EPSILON=$(echo ${NUM_NODES} | awk '{printf "%1.10f\n", 1 / ($1-1) }')
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
								### Skipping unnecessary iterations
								if [ "$NET_TYPE" != "erdos-renyi" ] && [ $LINK_PROBABILITY != ${LINK_PROBABILITY_LIST[0]} ];
								then
									continue
								fi
								if [ "$NET_TYPE" != "barabasi-albert" ] && [ $NUM_EDGES != ${NUM_EDGES_LIST[0]} ];
								then
									continue
								fi
								if [ "$UPDATE_CONF" == "belief-up" ] && [ $DEC_MODEL != "conf-perfect" ];
								then
									continue
								fi
								if [ $DEC_MODEL == "best-acc" ] && [ "$UPDATE_CONF" != "no-up" ];
								then
									continue
								fi
								
								#JOB_PARAM="net-${NET_TYPE}_nodes-${NUM_NODES}_edges-${NUM_EDGES}_link-${LINK_PROBABILITY}_model-${DEC_MODEL}_up-${UPDATE_CONF}_acc-${ACC_MEAN}_acstdv-${ACC_STD_DEV}_trunc-${ACC_TRUNCATED}"
								JOB_PARAM="net-${NET_TYPE}_nodes-${NUM_NODES}_edges-${NUM_EDGES}_link-${LINK_PROBABILITY}_model-${DEC_MODEL}_up-${UPDATE_CONF}_acc-${ACC_MEAN}_acstdv-${ACC_STD_DEV}"
								OUT_TXT="${OUTPUT_DATA_DIR}out_${JOB_PARAM}.txt"
								#OUT_PDF="${OUTPUT_DATA_DIR}pdf/out_${JOB_PARAM}"
									
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

echo "Submitted " $COUNT " jobs"
