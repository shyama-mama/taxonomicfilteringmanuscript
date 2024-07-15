#!/bin/bash


TOTAL_READS=20000000 # 20 Million

ENDO_FRACTIONS="0.1 1 2 4 10 20 30 40 50 60"


# Passing these as env variables
#HUMAN_FRACTION="0.1"
#BACT_FRACTION="0.9"
#REAGENT_FRACTION="0.1"

#ENDO=${prefix}.endogenous.simulated_reads_s1.fq.gz
#HUMAN=human.contamination.simulated_reads_s1.fq.gz
#BACT=bacterial.contamination.simulated_reads_s1.fq.gz
#REAGENT=reagent.contamination.simulated_reads_s1.fq.gz


for endo_fraction in $ENDO_FRACTIONS; do 

    saved_rand=$(echo $RANDOM)

    endo_reads=$(echo "scale=3; $TOTAL_READS * ($endo_fraction / 100) " | bc | cut -d\. -f1)
    remaining_reads=$(echo "scale=3; $TOTAL_READS-$endo_reads" | bc | cut -d\. -f1)
    human_reads=$(echo "scale=3; $remaining_reads * $HUMAN_FRACTION" | bc | cut -d\. -f1 )
    reagent_reads=$human_reads
    bact_reads=$(echo "scale=3; $TOTAL_READS-$endo_reads-$reagent_reads-$human_reads" | bc | cut -d\. -f1)  
    reagent_reads=$(echo "scale=3; $remaining_reads * $REAGENT_FRACTION" | bc | cut -d\. -f1)

    echo "endo_reads: ${endo_reads}"
    echo "human_reads: ${human_reads}"
    echo "bact_reads: ${bact_reads}"
    echo "reagent_reads: ${reagent_reads}"
    echo "random_seed: ${saved_rand}"

    seqtk sample -s${saved_rand} ${ENDO} ${endo_reads} > ${prefix}_seed${saved_rand}_endo${endo_fraction}.fq
    seqtk sample -s${saved_rand} ${HUMAN} ${human_reads} >> ${prefix}_seed${saved_rand}_endo${endo_fraction}.fq
    seqtk sample -s${saved_rand} ${BACT} ${bact_reads} >> ${prefix}_seed${saved_rand}_endo${endo_fraction}.fq
    seqtk sample -s${saved_rand} ${REAGENT} ${reagent_reads} >> ${prefix}_seed${saved_rand}_endo${endo_fraction}.fq

    gzip ${prefix}_seed${saved_rand}_endo${endo_fraction}.fq
done
