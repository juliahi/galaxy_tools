#!/bin/bash

USAGE="Usage: $0 annotation counts_output_file txt_output_file pdf_output_file proc_num [-e|-i] -b treated_plus treated_minus [t_plus2 t_minus2 ...] -u untreated_plus untreated_minus [...]"

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ANNOT=$1
OUTPUT1=$2
OUTPUT2=$3
OUTPUT3=$4
NPROC=$5
FITTYPE=$6
args=("$@")

echo `date -u`
#echo "Running with args:" ${args[@]}
OPTE=""
OPTI=""
OPTT=""

TREATEDP=()
UNTREATEDP=()
NAMEST=()
TREATEDM=()
UNTREATEDM=()
NAMESN=()
X=()

index=6
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
          NAME="${args[$index]}"
	    ((index++)) 
          TREATEDP=(${TREATEDP[@]} "${NAME}")
          IND="${args[$index]}"
	    ((index++)) 
          if [ "${IND}" != "" ] && [ ! "${NAME}.bai" ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi 
          NAME="${args[$index]}"
	    ((index++)) 
          TREATEDM=(${TREATEDM[@]} "${NAME}")
          IND="${args[$index]}"
	    ((index++)) 
          if [ "${IND}" != "" ] && [ ! "${NAME}.bai" ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi 
          NAMEST=(${NAMEST[@]} "'${args[$index]}'")
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
          NAME="${args[$index]}"
	    ((index++)) 
          UNTREATEDM=(${UNTREATEDM[@]} "${NAME}")
          IND="${args[$index]}"
	    ((index++)) 
          if [ ${IND} != "" ] && [ ! ${NAME}.bai ]; then
             ln -s "${IND}" "${NAME}.bai"
          fi 
          NAMESU=(${NAMESU[@]} "'${args[$index]}'")
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

#echo "treated: ${NAMEST} ${TREATEDP[@]} - ${TREATEDM[@]}"
#echo "untreated: ${NAMESU} ${UNTREATEDP[@]} - ${UNTREATEDM[@]}"


(python ${DIR}/count_for_deseq_par_list.py $ANNOT -s "+" $OPTE $OPTI $OPTT -o "${OUTPUT1}.plus" -b ${TREATEDP[@]} ${UNTREATEDP[@]} -n $(( NPROC / 2 )) -l ${NAMEST} ${NAMESU} ; echo $? >"$tmp1" ) &
proc1=$!

(python ${DIR}/count_for_deseq_par_list.py "$ANNOT" -s "-" $OPTE $OPTI $OPTT -o "${OUTPUT1}.minus" -b ${TREATEDM[@]} ${UNTREATEDM[@]} -n $(( NPROC / 2 )) -l ${NAMEST} ${NAMESU} ; echo $? >"$tmp2" ) &
proc2=$!
wait ${proc1} ${proc2}

read ret1 <"$tmp1"
read ret2 <"$tmp2"
rm "$tmp1" "$tmp2"

echo `date -u` "Finished counting reads"
if [ "$ret1" != 0 ] || [ "$ret2" != 0 ]; then
  >&2 echo "Counting reads failed. Exit status $ret1, $ret2"
  exit 1
fi

cat "${OUTPUT1}.plus" > "${OUTPUT1}"
wait
tail -n +2 "${OUTPUT1}.minus" >> "${OUTPUT1}"
wait


${DIR}/run_deseq.R "${OUTPUT1}" "${OUTPUT2}" "${OUTPUT3}" $FITTYPE ${X[@]} &
wait $!
x=$?


echo `date -u` "DESeq ended with code: $x"
exit $x

