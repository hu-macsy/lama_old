include ../make.inc

INCLUDE_FLAGS = -I${SCAI_BOOST_INCLUDE_DIR} -I${SCAI_ROOT}/include -I${SCAI_CUDA_INCLUDE_DIR}

LD_FLAGS = ${SCAI_LD_FLAGS} -Wl,-rpath=${SCAI_ROOT}/lib:${SCAI_CUDA_LIBRARY_PATH}

EXECUTABLE: ${EXECUTABLES_CUDA}

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .cu .o
.SUFFIXES: .o .exe

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS}
	
.cpp.o:
	${CXX} ${CXX_FLAGS} ${SCAI_DEFINES} ${INCLUDE_FLAGS} -c $*.cpp
	
.cu.o:
	${CXX_CUDA} ${CXX_CUDA_FLAGS} ${SCAI_DEFINES} ${INCLUDE_FLAGS} -c $*.cu

clean:
	rm -rf *.exe *.o
