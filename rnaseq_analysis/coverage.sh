#!/bin/bash

USAGE="Usage: $0 -o output -i bamfile index [bam2 id2 ...] [-n]\n"

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
args=("$@")

OPTN=""

BAMS=()

index=0
while [ $index -lt ${#args[@]} ]; do
  case ${args[$index]} in 
    "-n") OPTN="-n"
	    ((index++));; 
    "-o") 
	    ((index++)) 
      OUTPUT=${args[$index]}
	    ((index++));; 
    "-i") 
	    ((index++))
        while [ $index -lt ${#args[@]} ] && [ ${args[$index]:0:1} != "-" ]; do
          NAME=${args[$index]}
          BAMS=(${BAMS[@]} "${NAME}")
	    ((index++)) 
          IND=${args[$index]#?}
	    ((index++))
          if [ ${IND} != "" ] && [ ! -f "${NAME}.bai" ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi
       done;;
    *) echo "$USAGE" 1>&2
       echo "$index ${args[$index]}" 1>&2
       exit 1;;
  esac
done

# echo `date -u`

python ${DIR}/coverage.py ${OPTN} -o "${OUTPUT}" -i ${BAMS[@]}


if [ $? != 0 ] ; then
  echo "Coverage failed. Exit status $ret1"
  exit $ret1 
fi

