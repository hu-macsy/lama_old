include ../make.inc

INCLUDE_FLAGS = -I${BOOST_INCLUDE_DIR} -I${LAMA_ROOT}/include

LD_FLAGS = -L${LAMA_ROOT}/lib -lama -Wl,-rpath=${LAMA_ROOT}/lib 

EXECUTABLE: task0.exe task1a.exe task1b.exe task2.exe task2a.exe task3.exe task4.exe task5.exe

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS} 

.cpp.o:
	${CXX} ${CXX_FLAGS} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o *.frm *.frv *.amg *.vec *.trace
