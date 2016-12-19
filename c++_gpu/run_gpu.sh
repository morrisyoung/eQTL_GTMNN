#!/bin/bash
#$ -cwd -S /bin/bash -j y
##$ -pe smp 3
#$ -l mem=60G,time=160::
#$ -l gpu=1
#$ -M sy2515@c2b2.columbia.edu -m bes


./main 3 500 0.0000001
#./main_filefmt

