MPIF90=mpiifort 
CC=mpiicc -fPIC
FFLAGS1 = -traceback -O3 -fPIC

INC=-I./

MKLROOT=/opt/dario/intel/mkl

MKL_LIB=-L$(MKLROOT)/lib/intel64  -lmkl_blas95_lp64 -lmkl_lapack95_lp64  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread

SPARSKIT=libskit.a
CSPARSE=libcsparse.a
MCFIEDLER=libmcfiedler

# todo:
# sudo ln -s /usr/lib/libskit.so.2.0 /usr/lib/libskit.so
LD_FLAGS = -lskit -lcsparse $(MKL_LIB)

TEST_MATRIX=matrices/small.mtx
NUM_PROCS=4

###LIB_OBJS = sccblocks.o idxsort.o kb07ad.o preprocess.o laplacian.o matvec.o \

LIB_OBJS = sccblocks.o idxsort.o preprocess.o laplacian.o matvec.o \
fiedlerdense.o mc_general_fiedler.o

#========================================================================

all: mcfiedler_lib_shared


mcfiedler_lib: $(LIB_OBJS)
	ar r  $(MCFIEDLER).a $(LIB_OBJS)

mcfiedler_lib_shared: $(LIB_OBJS)
	$(MPIF90) -shared -o $(MCFIEDLER).so $(LIB_OBJS) $(LD_FLAGS)

test: $(LIB_OBJS) test.f90
	$(MPIF90) -o $@ test.f90 $(LIB_OBJS) $(LD_FLAGS)

tar:
	/bin/rm -r -f  TraceMIN_Fiedler
	mkdir TraceMIN_Fiedler
	cp Makefile TraceMIN_Fiedler
	cp cs.h TraceMIN_Fiedler
	cp fiedlerdense.f90 TraceMIN_Fiedler
	cp idxsort.c TraceMIN_Fiedler
	cp laplacian.f90  TraceMIN_Fiedler
	cp matvec.f90  TraceMIN_Fiedler
	cp mc_general_fiedler.f90  TraceMIN_Fiedler
	cp preprocess.f90  TraceMIN_Fiedler
	cp pymcfiedlermodule.c  TraceMIN_Fiedler
	cp sccblocks.c TraceMIN_Fiedler
	tar cvfz TraceMIN_Fiedler_050313.tar.gz TraceMIN_Fiedler
	/bin/rm -r TraceMIN_Fiedler/*
	rmdir TraceMIN_Fiedler

clean:
	/bin/rm -f *.exe *.o *.so

#========================================================================

.SUFFIXES:
.SUFFIXES: .f90 .f .c .o

.f90.o:
	$(MPIF90) -c  $(FFLAGS1) $(INC)  $< 
.f.o:
	$(MPIF90) -c  $(FFLAGS1) $(INC)  $< 
.c.o:
	$(CC) -c $(INC)  $<

#========================================================================
