include ../make.inc

INCLUDE_FLAGS = -I${BOOST_INCLUDE_DIR} -I${LAMA_ROOT}/include

LD_FLAGS = ${LAMA_LD_FLAGS} -Wl,-rpath=${LAMA_ROOT}/lib 

EXECUTABLE: simple.exe vector.exe solver.exe blas1.exe vector_exp.exe matrix.exe scalar.exe spmvCPUCUDA.exe cgCPUCUDA.exe

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS} 

.cpp.o:
	${CXX} ${CXX_FLAGS} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o *.frm *.frv *.amg *.vec *.trace
