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
NORMS=()

index=0
while [ $index -lt ${#args[@]} ]; do
  case ${args[$index]} in 
    "-o") 
            ((index++)) 
      OUTPUT1=${args[$index]}
            ((index++)) 
      OUTPUT2=${args[$index]}
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
          
          NAME=${args[$index]}
          BAMS=(${BAMS[@]} "${NAME}")
            ((index++)) 
          IND=${args[$index]#?}
            ((index++))
          if [ ${IND} != "" ] && [ ! -f "${NAME}.bai" ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi
          
          TYPE=${args[$index]}
            ((index++))
          if [ ${TYPE} == "no" ]; then
             NORMS=(${NORMS[@]} 1000000)
          elif [ ${TYPE} == "reads" ]; then
             NORMS=(${NORMS[@]} "`samtools view -c ${NAME}`")
          elif [ ${TYPE} == "number" ]; then
             NORMS=(${NORMS[@]} "${args[$index]}")
            ((index++))
          elif [ ${TYPE} == "file" ]; then
             NORMS=(${NORMS[@]} "`cat ${args[$index]}`")
            ((index++))
          else 
             echo "no type for file $NAME"
             exit 1
          fi
       done;;
    *) echo "$USAGE" 1>&2
       echo "$index ${args[$index]}" 1>&2
       exit 1;;
  esac
done

# echo `date -u`

python ${DIR}/coverage_paired.py  -o "${OUTPUT1}" "${OUTPUT2}"  -i "${BAMS[@]}" -n "${NORMS[@]}"
echo ${DIR}/coverage_paired.py  -o "${OUTPUT1}" "${OUTPUT2}"  -i "${BAMS[@]}" -n "${NORMS[@]}"

if [ $? != 0 ] ; then
  echo "Coverage failed. Exit status $ret1"
  exit $ret1 
fi

