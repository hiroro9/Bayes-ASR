#!/bin/bash
#============ LSF Options ============
#QSUB -q gr10083b
#QSUB -ug gr10083
#QSUB -W 300:00
#QSUB -A p=1:t=3:c=3:m=10G
#QSUB -o output/model-3-data-01-15000.out
#QSUB -e output/model-3-data-01-15000.err
#============ Shell Script ============
set -x

file="model-3-data-01-15000"

#=========================================#
log=".log"
output="./output/"

~/R-3.6.3-sysB/bin/Rscript model-3.R > $output$file$log
