#!/bin/bash
set -e

make clean
rm -rf lama.png
rm -rf gr_30_30.mtx.*.png
rm -rf *.exe *.o *.frm *.mtx *.frv *.amg *.vec *.trace

# build examples
make

# run examples

# lecture/*
lecture/task0.exe lecture/gr_30_30.mtx
lecture/task1a.exe lecture/gr_30_30.mtx
lecture/task1b.exe
lecture/task2.exe lecture/gr_30_30.mtx
lecture/task2a.exe lecture/gr_30_30.mtx
lecture/task3.exe lecture/gr_30_30.mtx
lecture/task4.exe lecture/gr_30_30.mtx
lecture/task5.exe lecture/gr_30_30.mtx

# check if there are unkown examples
count=`ls -l -la lecture/*.exe | wc -l`
if [ $count -ne 8 ]; then
	echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi

# solver/*
cd solver
./matrix_generator.exe example 3 27 100 100 100
./matrix_convert.exe -mm example.frv example.mtx
./vector_generator.exe example2.mtx 1000 1
./cg_solver.exe example
./gmres_solver.exe example
#./jacobi_solver.exe example
#./amg_solver.exe example
./lama_info.exe
cd ..

# check if there are unkown examples
count=`ls -l -la solver/*.exe | wc -l`
if [ $count -ne 8 ]; then
	echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi

# spy/*
spy/amg_spy.exe lecture/gr_30_30.mtx
spy/spy.exe lecture/gr_30_30.mtx

# check if there are unkown examples
count=`ls -l -la spy/*.exe | wc -l`
if [ $count -ne 2 ]; then
	echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi
