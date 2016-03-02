include make.inc

INCLUDE_FLAGS = -I${SCAI_BOOST_INCLUDE_DIR} -I${SCAI_ROOT}/include

LD_FLAGS = ${SCAI_LD_FLAGS} -Wl,-rpath=${SCAI_ROOT}/lib 

EXECUTABLE: ${EXECUTABLES} Module.so DummyModule.so

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe
.SUFFIXES: .o .so

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS}

.o.so:
	${CXX} -o lib$*.so $*.o ${LD_FLAGS} -shared
	
.cpp.o:
	${CXX} -fPIC ${CXX_FLAGS} ${SCAI_DEFINES} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o *.so