#!/bin/bash

USAGE="Usage: $0 output s window step slots [-e|-i] -b bamfile index [bam2 id2 ...] "

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

OUTPUT1=$1
S=$2
WINDOW=$3
STEP=$4
SLOTS=$5
args=("$@")

echo ${args[@]}

BAMS=()
NAMES=()

index=5
while [ $index -lt ${#args[@]} ]; do
  case ${args[$index]} in 
    "-b") 
	    ((index++))
        while [ $index -lt ${#args[@]} ] && [ ${args[$index]:0:1} != "-" ]; do
          NAME=${args[$index]}
          BAMS=(${BAMS[@]} "${NAME}")
	    ((index++)) 
          IND=${args[$index]#?}
	    ((index++))
          NAMES=(${NAMES[@]} "${args[$index]}")
	    ((index++))
          if [ ${IND} != "" ] && [ ! ${NAME}.bai ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi
          #echo "${NAME} ${IND}"
       done;;
    *) echo "$USAGE"
       echo "$index ${args[$index]}"
       exit 1;;
  esac
done

tmp1=$(mktemp)


echo `date -u`

(python ${DIR}/count_window.py -s $S --window ${WINDOW} --step ${STEP} -o "${OUTPUT1}" -b ${BAMS[@]} -l ${NAMES[@]} -n ${SLOTS} ; echo $? >"$tmp1" ) &
proc1=$!

wait ${proc1} 

read ret1 <"$tmp1"
rm "$tmp1" 


echo "Finished: " `date -u`

if [ "$ret1" != 0 ] ; then
  echo "Counting reads failed. Exit status $ret1"
  exit $ret1 
fi

exit $x

