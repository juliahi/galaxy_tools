#!/bin/bash

USAGE="Usage: $0 annotation output_file proc_num overlap [+/-]  -f [f f.bai label ...] "

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ANNOT=$1
OUTPUT1=$2
NPROC=$3
OVERLAP=$4
STRAND=$5
args=("$@")

echo `date -u`
#echo "Running with args:" "$@"

TREATEDP=()
UNTREATEDP=()
NAMEST=()

index=5

while [ $index -lt ${#args[@]} ]; do
  case "${args[$index]}" in 
    "-f") 
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
        done;;
    *) >&2 echo "$USAGE"
       >&2 echo "$index ${args[$index]}"
       exit 1;;
  esac
done

tmp1=$(mktemp)
tmp2=$(mktemp)

echo files: "${NAMEST[@]}" ${TREATEDP[@]} - ${TREATEDM[@]}
#echo untreated: ${NAMESU[@]} ${UNTREATEDP[@]} - ${UNTREATEDM[@]}

python ${DIR}/policystr_annotRom_mate.py -a "$ANNOT" -s "${STRAND}" -o "${OUTPUT1}.plus" -f ${TREATEDP[@]}  -c "${NAMEST[@]}"  -n $OVERLAP; echo $? >"$tmp1" 


read ret1 <"$tmp1"
rm "$tmp1" 

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


