all: makeAll

clean: makeAllClean

DIRS = bench \
	   labelrank \
	   lecture \
       solver \
       spy \
       tutorial

makeAll:
	for dir in $(DIRS); do \
		( cd $$dir; make ); \
	done;

makeAllClean:
	for dir in $(DIRS); do \
		( cd $$dir; make clean ); \
	done;

