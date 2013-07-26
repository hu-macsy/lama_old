all: makeAll

clean: makeAllClean

DIRS = lecture \
       tutorial \
       bench \
       biCGTest \
       solver

makeAll:
	for dir in $(DIRS); do \
		( cd $$dir; make ); \
	done;

makeAllClean:
	for dir in $(DIRS); do \
		( cd $$dir; make clean ); \
	done;

