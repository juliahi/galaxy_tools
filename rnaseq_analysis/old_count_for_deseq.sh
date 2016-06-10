#!/bin/bash

USAGE="Usage: $0 annotation output s slots [-e|-i] -b bamfile index [bam2 id2 ...] "

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ANNOT=$1
OUTPUT1=$2
S=$3
SLOTS=$4
args=("$@")

echo "${args[@]}"
OPTE=""
OPTI=""
OPTT=""

BAMS=()
NAMES=()

index=4
while [ $index -lt ${#args[@]} ]; do
  case ${args[$index]} in 
    "-e") OPTE="-e"
	    ((index++));; 
    "-i") OPTI="-i"
	    ((index++));; 
    "-t") OPTT="-t"
	    ((index++));; 
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
if [ "$OPTE" == "-e" ]; then
  echo "Starting count of exons"
else 
  if [ "$OPTT" == "-t" ]; then
    echo "Starting count of transcripts"
  else
    echo "Starting count of genes"
  fi  
fi

#(python ${DIR}/count_for_deseq_par_list.py $ANNOT -s $S $OPTE $OPTI $OPTT -o "${OUTPUT1}" -b "${BAMS[@]}" -l "${NAMES[@]}" -n ${SLOTS} ; echo $? >"$tmp1" ) &
#proc1=$!
echo running python ${DIR}/count_for_deseq_par_list.py $ANNOT -s $S $OPTE $OPTI $OPTT -o "${OUTPUT1}" -b "${BAMS[@]}" -l "${NAMES[@]}" -n ${SLOTS}
exit
wait ${proc1} 

read ret1 <"$tmp1"
rm "$tmp1" 


echo "Finished: " `date -u`

if [ "$ret1" != 0 ] ; then
  echo "Counting reads failed. Exit status $ret1"
  exit $ret1 
fi

exit $x

