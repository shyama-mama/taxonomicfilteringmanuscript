#!/bin/bash

# Size distribution and matrix file
FREQUENCY_FILE=gargammel/src/sizefreq.size.gz
MATRIX_PREFIX=gargammel/src/matrices/single-

# Simulate modern human contamination
gargammel -n 20000000  --comp 0,1,0 \
    -f ${FREQUENCY_FILE} \
    -ss HSXt \
    -o human.contamination.simulated_reads human/

# Simulate reagent contamination
gargammel -n 20000000  --comp 1,0,0 \
    -f ${FREQUENCY_FILE} \
    -ss HSXt \
    -o reagent.contamination.simulated_reads reagent/

# simulate ancient bacteria
gargammel -n 20000000  --comp 1,0,0 \
    -f ${FREQUENCY_FILE} \
    -matfile ${MATRIX_PREFIX} \
    -ss HSXt \
    -o bacterial.contamination.simulated_reads bacteria/

# simulate endogenous human and dog
for species in $(echo "human dog"); do
    gargammel -n 20000000  --comp 1,0,0 \
        -f ${FREQUENCY_FILE} \
        -matfile ${MATRIX_PREFIX} \
        -ss HSXt \
        -o ${species}.endogenous.simulated_reads ${species}/
done
