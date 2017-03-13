#!/bin/bash

S1='string'
S2='string'
S3=$S1

if [ "$S1" == "$S2" ];
then
        echo "S1('$S1') is equal to S2('$S2')"
fi
if [ "$S1" == "$S3" ];
then
        echo "S1('$S1') is equal to S3('$S3')"
fi

NET_TYPE_LIST=('full' 'erdos-renyi' 'barabasi-albert')
NUM_EDGES_LIST=($(seq 2 3 9))

#echo ${NUM_EDGES_LIST[*]}
#echo ${NUM_EDGES_LIST[0]}
#echo ${NET_TYPE_LIST[0]}

for NET_TYPE in ${NET_TYPE_LIST[*]}
do
	for NUM_EDGES in ${NUM_EDGES_LIST[*]}
	do
		echo $NET_TYPE
		echo $NUM_EDGES
		if [ "$NET_TYPE" != "barabasi-albert" ] && [ $NUM_EDGES != ${NUM_EDGES_LIST[0]} ];
		then
			echo "Got not a B-A and here I should skip"
		fi
	done
done
