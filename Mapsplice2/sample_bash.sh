#! /bin/bash

OUTPUT_DIR=/home/small_results/ 
READ_FILE_END1=/home/small_1.fastq
READ_FILE_END2=/home/small_2.fastq
MAPSPLICE_DIR=/home/MapSplice-v2.1.3/
REF_GENOME=/home/hg19_sequence/
BOWTIE_INDEX=/home/hg19_index/humanchridx_M


python $MAPSPLICE_DIR/mapsplice.py \
       -1 $READ_FILE_END1 \
       -2 $READ_FILE_END2 \
       -c $REF_GENOME \
       -x $BOWTIE_INDEX \
       -p 10 \
       -o $OUTPUT_DIR 2>log.txt
 

