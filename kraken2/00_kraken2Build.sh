#!/bin/bash

# Assumes you've already download the needed libraries to build the database. 

# defaults (k-mer 35):
#kmer_len=35
#minimizer_len=32
#minimizer_spaces=7

# aDNA optamised (k-mer 29):
#kmer_len=29
#minimizer_len=24
#minimizer_spaces=6

kraken2-build --build --db $DBNAME --threads 64 --kmer-len ${kmer_len} --minimizer-len ${minimizer_len} --minimizer-spaces ${minimizer_spaces}
