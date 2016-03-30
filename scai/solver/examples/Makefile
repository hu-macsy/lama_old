all: makeAll

clean: makeAllClean

DIRS = lecture \
	   myJacobi \
	   solver \
       spy

makeAll:
	for dir in $(DIRS); do \
		( cd $$dir; make ); \
	done;

makeAllClean:
	for dir in $(DIRS); do \
		( cd $$dir; make clean ); \
	done;

