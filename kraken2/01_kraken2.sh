#!/bin/bash

kraken2 \
    --db ${DB} \
    --output ${prefix}.kraken.out \
    --report ${prefix}.kraken.report \
    --confidence ${CONFIDENCE} \
    --gzip-compressed ${fastq}
