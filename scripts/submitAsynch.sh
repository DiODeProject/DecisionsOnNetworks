#!/bin/bash

PROJECT_HOME="${HOME}/DecisionsOnNetworks"
CONF_DIR="${PROJECT_HOME}/conf_cluster"
mkdir -p ${CONF_DIR}
TEMPLATE_SETTINGS="${PROJECT_HOME}/conf/AsynchK.template.config"
PYPATH="${PROJECT_HOME}/src/"
EXEC_FILE="${PROJECT_HOME}/src/AsynchKicks/DecisionProcess.py"
OUTPUT_DATA_DIR="${PROJECT_HOME}/kick_5/"

SEED=221189
NUM_EXP=100
NUM_RUN=1
MAX_TIME=1000

AGENT_TYPE='DDM'
UPDATE_MODEL_LIST=( 'no-up' 'thresh-kick' 'conf-kick' )

### Unnecessary for this analysis 
DDM_DRIFT_DIST='normal'
DDM_BASE_DRIFT_LIST=( 0.05 0.1 0.2 )
#DDM_DRIFT_RANGE_LIST=( 0.25 0.50 1.00 )
DDM_DRIFT_RANGE=0.5
#DDM_DRIFT_STD_DEV_LIST=( 0.25 0.50 1.00 )
#DDM_DRIFT_STD_DEV_LIST=( 0.10 0.15 0.20 0.30 0.40 )
DDM_DRIFT_STD_DEV_LIST=($(seq 0.10 0.10 0.50))
DDM_NOISE_STD_DEV=0.5
THRESHOLD=10
PRIOR_DIST=0.5
USE_BAYES_RISK='true'
COST_MATRIX_T=1
COST_MATRIX_E=20

#NUM_NODES=20
NUM_NODES_LIST=($(seq 20 10 100))
NET_TYPE_LIST=('erdos-renyi' 'barabasi-albert' 'space')
#NET_TYPE_LIST=('space')
#NET_TYPE_LIST=('erdos-renyi' 'barabasi-albert')
#LINK_PROBABILITY_LIST=($(seq 0.2 0.1 0.8))
#LINK_PROBABILITY_LIST=(0.2 0.4 0.6 0.8)
LINK_PROBABILITY_LIST=(0.2 0.6)
#NUM_EDGES=3
#NUM_EDGES_LIST=($(seq 3 3 12))
NUM_EDGES_LIST=(3 9)
#INTERACTION_RADIUS_LIST=($(seq 0.05 0.05 0.3))
#INTERACTION_RADIUS_LIST=($(seq 0.20 0.05 0.35))
INTERACTION_RADIUS_LIST=( 0.20 )
ENV_SIZE=1
#ENV_SIZE_LIST=(0.5 2 3 4 5 6 7 8 9 15 20)
#ENV_SIZE_LIST=(0.01 0.1 0.5 1 2 3 4 5 10 15 20)

DYNAMICNET='false'
#AGENTSPEED_LIST=($(seq 0.05 0.05 0.2))
AGENTSPEED=0.05 
MOVE_STDDEV=45
PERIODICBOUND='false'

ACC_MEAN=0.6
ACC_STD_DEV=0.12
ACC_TRUNCATED='true'

COUNT=0

for UPDATE_MODEL in ${UPDATE_MODEL_LIST[*]}
do
	for NET_TYPE in ${NET_TYPE_LIST[*]}
	do
		for NUM_NODES in ${NUM_NODES_LIST[*]}
		do
			for INTERACTION_RADIUS in ${INTERACTION_RADIUS_LIST[*]}
			do
				for LINK_PROBABILITY in ${LINK_PROBABILITY_LIST[*]}
				do
					for NUM_EDGES in ${NUM_EDGES_LIST[*]}
					do
						for DDM_DRIFT_STD_DEV in ${DDM_DRIFT_STD_DEV_LIST[*]}
						do
							for DDM_BASE_DRIFT in ${DDM_BASE_DRIFT_LIST[*]}
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
								if [ "$NET_TYPE" != "space" ] && [ $INTERACTION_RADIUS != ${INTERACTION_RADIUS_LIST[0]} ];
								then
									continue
								fi
								
								NET_PAR=0
								if [ "$NET_TYPE" == "erdos-renyi" ];
								then
									NET_PAR=${LINK_PROBABILITY}
								fi
								if [ "$NET_TYPE" == "barabasi-albert" ];
								then
									NET_PAR=${NUM_EDGES}
								fi
								if [ "$NET_TYPE" == "space" ];
								then
									NET_PAR=${INTERACTION_RADIUS}
								fi
								
								JOB_PARAM="net-${NET_TYPE}_nodes-${NUM_NODES}_link-${NET_PAR}_up-${UPDATE_MODEL}_driftbase-${DDM_BASE_DRIFT}_driftrange-${DDM_DRIFT_STD_DEV}"
								OUT_TXT="${OUTPUT_DATA_DIR}out_${JOB_PARAM}.txt"
									
								CONF_FILE="${CONF_DIR}/kick_${JOB_PARAM}.config"
								
								sed -e "s|SEED|${SEED}|" \
									-e "s|NUM_EXP|${NUM_EXP}|" \
									-e "s|NUM_RUN|${NUM_RUN}|" \
									-e "s|MAX_TIME|${MAX_TIME}|" \
									-e "s|AGENT_TYPE|${AGENT_TYPE}|" \
									-e "s|ACC_MEAN|${ACC_MEAN}|" \
									-e "s|ACC_STD_DEV|${ACC_STD_DEV}|" \
									-e "s|ACC_TRUNCATED|${ACC_TRUNCATED}|" \
									-e "s|UPDATE_MODEL|${UPDATE_MODEL}|" \
									-e "s|DDM_DRIFT_DIST|${DDM_DRIFT_DIST}|" \
									-e "s|DDM_BASE_DRIFT|${DDM_BASE_DRIFT}|" \
									-e "s|DDM_DRIFT_RANGE_PLUS_MINUS|${DDM_DRIFT_RANGE}|" \
									-e "s|DDM_DRIFT_STD_DEV|${DDM_DRIFT_STD_DEV}|" \
									-e "s|DDM_NOISE_STD_DEV|${DDM_NOISE_STD_DEV}|" \
									-e "s|THRESHOLD|${THRESHOLD}|" \
									-e "s|PRIOR_DIST|${PRIOR_DIST}|" \
									-e "s|USE_BAYES_RISK|${USE_BAYES_RISK}|" \
									-e "s|COST_MATRIX_T|${COST_MATRIX_T}|" \
									-e "s|COST_MATRIX_E|${COST_MATRIX_E}|" \
									-e "s|NET_TYPE|${NET_TYPE}|" \
									-e "s|NUM_NODES|${NUM_NODES}|" \
									-e "s|NUM_EDGES|${NUM_EDGES}|" \
									-e "s|LINK_PROBABILITY|${LINK_PROBABILITY}|" \
									-e "s|INTERACTION_RADIUS|${INTERACTION_RADIUS}|" \
									-e "s|DYNAMICNET|${DYNAMICNET}|" \
									-e "s|AGENTSPEED|${AGENTSPEED}|" \
									-e "s|MOVE_STDDEV|${MOVE_STDDEV}|" \
									-e "s|ENV_SIZE|${ENV_SIZE}|" \
									-e "s|PERIODICBOUND|${PERIODICBOUND}|" \
									-e "s|OUT_TXT|${OUT_TXT}|" \
										${TEMPLATE_SETTINGS} > ${CONF_FILE}
											
								export PYTHONPATH=${PYPATH}
								#COMMAND="python3 ${EXEC_FILE} ${CONF_FILE}"
								#COMMAND="./run_job.sh ${EXEC_FILE} ${CONF_FILE}"
								COMMAND="qsub run_job.sh ${EXEC_FILE} ${CONF_FILE}"
								${COMMAND}
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
