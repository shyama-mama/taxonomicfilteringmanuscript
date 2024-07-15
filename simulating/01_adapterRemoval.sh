#!/bin/bash

set -e 
set -x 

mkdir -p ${OUTPUT}

AdapterRemoval --file1 $R1 --file2 $R2 \
    --basename $PREFIX --gzip --threads 2 \
    --qualitymax 64 --collapse  --trimns --trimqualities \
    --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
    --minlength 30 --minquality 20 --minadapteroverlap 1

cat ${PREFIX}.collapsed.gz ${PREFIX}.collapsed.truncated.gz ${PREFIX}.singleton.truncated.gz ${PREFIX}.pair1.truncated.gz ${PREFIX}.pair2.truncated.gz > ${OUTPUT}/${PREFIX}.combined.tmp.fq.gz

mv ${PREFIX}.settings ${OUTPUT}/

## Add R_ and L_ for unmerged reads for DeDup compatibility
AdapterRemovalFixPrefix -Xmx4g ${OUTPUT}/${PREFIX}.combined.tmp.fq.gz | pigz -p 1 > ${OUTPUT}/${PREFIX}.combined.fq.gz
