all: makeAll

clean: makeAllClean

DIRS = bench \
	   labelrank \
       tutorial

makeAll:
	for dir in $(DIRS); do \
		( cd $$dir; make ); \
	done;

makeAllClean:
	for dir in $(DIRS); do \
		( cd $$dir; make clean ); \
	done;

