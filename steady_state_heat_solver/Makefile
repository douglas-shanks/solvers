
F95       = mpif90 
FLAGS     =
LIB_LIST  =-llapack -lblas
all: test_steady

greetings.o: greetings.f90
	$(F95) $(OPT) $(FLAGS) -c greetings.f90

header.o: header.f90
	$(F95) $(OPT) -c header.f90

sparsegather.o: sparsegather.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c sparsegather.f90

sysmatrix.o: sysmatrix.f90 header.o
	$(F95) $(OPT) -c sysmatrix.f90

sysmatrix_pre.o: sysmatrix_pre.f90 header.o
	$(F95) $(OPT) -c sysmatrix_pre.f90

rhs.o: rhs.f90 header.o
	$(F95) $(OPT) -c rhs.f90
	
cg.o: cg.f90 header.o
	$(F95) $(OPT) -c cg.f90
	
alpha.o: alpha.f90 header.o
	$(F95) $(OPT) -c alpha.f90	

vecdot.o: vecdot.f90 header.o
	$(F95) $(OPT) -c vecdot.f90	

pcg.o: pcg.f90 header.o
	$(F95) $(OPT) -c pcg.f90	
	
error.o: error.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c error.f90

matmult.o: matmult.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c matmult.f90

jacobi.o: jacobi.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c jacobi.f90

FDsolution.o: FDsolution.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c FDsolution.f90

test_steady.o: test_steady.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c test_steady.f90

test_steady: greetings.o sparsegather.o header.o error.o matmult.o sysmatrix_pre.o sysmatrix.o rhs.o cg.o pcg.o jacobi.o test_steady.o alpha.o vecdot.o FDsolution.o
	 $(F95) $(OPT) -o test_steady greetings.o header.o sparsegather.o error.o matmult.o sysmatrix.o sysmatrix_pre.o rhs.o cg.o pcg.o jacobi.o test_steady.o alpha.o vecdot.o FDsolution.o $(LIB_LIST)

clean: 
	rm -rf *~ *.o *.mod core  
