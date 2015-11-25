include make.inc

INCLUDE_FLAGS = -I${BOOST_INCLUDE_DIR} -I${SCAI_ROOT}/include

LD_FLAGS = ${SCAI_LD_FLAGS} -Wl,-rpath=${SCAI_ROOT}/lib 

EXECUTABLE: Barrier.exe CriticalRegion.exe DemoFactory.exe DemoFactory1.exe DemoPointer.exe ExceptionDemo.exe UseModule.exe

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe
.SUFFIXES: .o .so

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS} 
	
.cpp.o:
	${CXX} ${CXX_FLAGS} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o *.so