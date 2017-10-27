all: makeAll EXECUTABLE

clean: makeAllClean

DIRS = cuda

makeAll:
	for dir in $(DIRS); do \
		( if [ -d $$dir ]; then ( cd $$dir; make; ); fi ) \
	done;

makeAllClean:
	for dir in $(DIRS); do \
		( if [ -d $$dir ]; then ( cd $$dir; make clean; ); fi ) \
	done;


include make.inc

INCLUDE_FLAGS = -I${SCAI_BOOST_INCLUDE_DIR} -I${SCAI_ROOT}/include

LD_FLAGS = ${SCAI_LD_FLAGS} -Wl,-rpath=${SCAI_ROOT}/lib 

EXECUTABLE: ${SCAI_EXECUTABLES}

.SUFFIXES:

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe

.o.exe: 
	${CXX} -o $*.exe $*.o ${LD_FLAGS}
	
.cpp.o:
	${CXX} -fPIC ${CXX_FLAGS} ${SCAI_DEFINES} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o
