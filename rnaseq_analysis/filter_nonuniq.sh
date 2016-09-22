#!/bin/bash

USAGE="Usage: $0 inputfile input_index output1 output2 "

if [ "$#" == "0" ]; then
        echo "$USAGE"
        exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

IND=$2
NAME=$1
if [ "${IND}" != "" ] && [ ! "${NAME}.bai" ]; then
             ln -s "${IND}" "${NAME}.bai"
fi 


tmp1=$(mktemp)

python ${DIR}/filter_bam_NonUniq_dupl.py -f $NAME -o $3 $4; echo $? >"$tmp1" 


read ret1 <"$tmp1"
rm "$tmp1" 

if [ "$ret1" != 0 ] ; then
  >&2 echo "Error filtering reads. Exit status $ret1"
  exit 1
fi
