#!/bin/bash
# submit me with "qsub /home/jiri/workspace/LAMA/src/framework/LAMAMV/JirisBenchFiles/bull_bench_job.sh"
#Run in cuda queue
#PBS -q cuda
#Allocate 8 nodes
#PBS -lnodes=8
#Start after 8 pm
#PBS -a 2000
#Join output an error stream
#PBS -j oe
#Write Output to file
#PBS -o /home/jiri/tmp/PCG_bull_scaling/bull_bench.job.pbs.log

cd /home/jiri/workspace/LAMA/build/bull/gcc/openmpi/Release/framework

./BenchmarkRunner --output=/home/jiri/tmp/PCG_bull_scaling/bull_PCG_scaling_Distributed.csv -f /home/jiri/workspace/LAMA/src/framework/LAMAMV/JirisBenchFiles/LAMAPCGBenchs_bull_scaling.beru &> /home/jiri/tmp/PCG_bull_scaling/bull_bench.job2.log