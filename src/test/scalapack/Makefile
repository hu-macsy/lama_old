

MKL_DIR=/home/brandes/intel/mkl/lib/intel64
MPI_DIR=/home/brandes/local/openmpi-1.4.3/lib

OBJS = pdinvdriver.o pdinvinfo.o pmatgen.o pmatgeninc.o pdinvchk.o

pdinvdriver: $(OBJS)
	ifort -openmp $(OBJS) -L$(MKL_DIR) -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -L$(MPI_DIR) -lmpi

pdinvdriver.o: pdinvdriver.f
	ifort -c pdinvdriver.f

pdinvinfo.o: pdinvinfo.f
	ifort -c pdinvinfo.f

pdinvchk.o: pdinvchk.f
	ifort -c pdinvchk.f

pmatgen.o: pmatgen.f
	ifort -c pmatgen.f

pmatgeninc.o: pmatgeninc.f
	ifort -c pmatgeninc.f

