include ../make.inc

INCLUDE_FLAGS = -I${BOOST_INCLUDE_DIR} -I${LAMA_ROOT}/include

LD_FLAGS = -L${LAMA_ROOT}/lib -llama -Wl,-rpath=${LAMA_ROOT}/lib 

EXECUTABLE: matrix_generator.exe cg_solver.exe jacobi_solver.exe

cg_solver.exe: LAMAConfig.hpp
jacobi_solver.exe: LAMAConfig.hpp

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS} 

.cpp.o:
	${CXX} ${CXX_FLAGS} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o *.frm *.frv *.amg *.vec *.trace
