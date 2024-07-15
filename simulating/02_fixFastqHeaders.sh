#!/bin/bash 

addNumber() {
    header=$1
    seq=$2
    quality=$3
    n=$4
    number=$(printf "%09d" $n)
    echo "${header}" | sed 's/$/:'$number'/' | awk -F  ':' '{print $1":"$2":"$3":"$4":"$6":"$5;}'
    echo $seq
    echo "+"
    echo $quality
}

mkdir -p ${OUTPUT_DIR}
rm -f ${OUTPUT_DIR}/tmp*fq

N=20
n=0
(
zcat $FASTQ | while read header && read sequence && read plus && read quality ; do
    ((i=i%N)); ((i++==0)) && wait
    addNumber $header $sequence $quality $n >> ${OUTPUT_DIR}/tmp.$i.fq & 
    n=$(($n+1))
done
)

cat ${OUTPUT_DIR}/tmp*fq | gzip - > ${OUTPUT}
