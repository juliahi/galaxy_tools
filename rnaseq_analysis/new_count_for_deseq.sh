#!/bin/bash

USAGE="Usage: $0 annotation output_file proc_num strand [-e|-i] -b treated1.bam treated1.bai name1 [...] -u untreated1.bam untreated1.bai name1 [...]"

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ANNOT=$1
OUTPUT1=$2
NPROC=$3
STRAND=$4
args=("$@")

echo `date -u`
#echo "Running with args:" "$@"
OPTE=""
OPTI=""
OPTT=""

TREATEDP=()
UNTREATEDP=()
NAMEST=()
NAMESN=()
X=()

index=4


while [ $index -lt ${#args[@]} ]; do
  case "${args[$index]}" in 
    "-f")  ((index++)) 
	OPTF="-f ${args[$index]}"
	((index++));; 
    "-i")  ((index++)) 
	OPTI="-i ${args[$index]}"
	((index++));; 
    "-b") 
	    ((index++))
        while [ $index -lt ${#args[@]} ] && [ ${args[$index]:0:1} != "-" ]; do
          NAME="${args[$index]}"
	    ((index++)) 
          TREATEDP=(${TREATEDP[@]} "${NAME}")
          IND="${args[$index]}"
	    ((index++)) 
          if [ "${IND}" != "" ] && [ ! "${NAME}.bai" ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi 
          NAMEST=("${NAMEST[@]}" "${args[$index]}")
	    ((index++))
          X=(${X[@]} 1)
        done;;
    "-u") 
	    ((index++)) 
        while [ $index -lt ${#args[@]} ] && [ ${args[$index]:0:1} != "-" ]; do
          NAME="${args[$index]}"
	    ((index++)) 
          UNTREATEDP=(${UNTREATEDP[@]} "${NAME}")
          IND="${args[$index]}"
	    ((index++)) 
          if [ ${IND} != "" ] && [ ! ${NAME}.bai ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi 
          NAMESU=("${NAMESU[@]}" "${args[$index]}")
	    ((index++)) 
          X=(${X[@]} 0)
        done;;
    *) >&2 echo "$USAGE"
       >&2 echo "$index ${args[$index]}"
       exit 1;;
  esac
done

tmp1=$(mktemp)
tmp2=$(mktemp)

#echo treated: "${NAMEST[@]}" ${TREATEDP[@]} - ${TREATEDM[@]}
#echo untreated: ${NAMESU[@]} ${UNTREATEDP[@]} - ${UNTREATEDM[@]}

python ${DIR}/count_for_deseq_par_list.py $ANNOT -s "$STRAND" $OPTF $OPTI -o "${OUTPUT1}.plus" -b ${TREATEDP[@]} ${UNTREATEDP[@]} -n $NPROC -l "${NAMEST[@]}" "${NAMESU[@]}" ; echo $? >"$tmp1"  
read ret1 <"$tmp1"

echo `date -u` "Finished counting reads"
if [ "$ret1" != 0 ] ; then
  >&2 echo "Counting reads failed. Exit status $ret1"
  exit 1
fi

(head -n 1 "${OUTPUT1}.plus" && tail -n +2 "${OUTPUT1}.plus" | sort -u) > "${OUTPUT1}"
wait

grep "SUPV" "${OUTPUT1}"
wait

rm "${OUTPUT1}.plus" 
wait


