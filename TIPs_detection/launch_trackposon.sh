#!/bin/bash

list_file=$(ls /home/database/Beech_genomes/Massane/RAWDATA/*1_H2VCTDSXY*)

for file in $list_file
do
qsub -v file=$file /home/abirami/MASSANE/TRACKPOSON.sh
done

