LIB_LIST = -llapack -lblas

# Above is the location of the blas library which is called by lublas2
# The filename with full path for this library is /usr/lib/libblas.a
# Since it is in /usr/lib we do not need to specify a path.
# The compiler looks there automatically.
#F95 = gfortran
F95 = mpif90
OPT = -O

all:test2D

header.o: header.f90
	$(F95) $(OPT) -c header.f90
	
sparsegather.o: sparsegather.f90 header.o
	$(F95) $(OPT) -c sparsegather.f90
	
initial2D.o: initial2D.f90 header.o
	$(F95) $(OPT) -c initial2D.f90	
	
Mat_Mult.o: Mat_Mult.f90 header.o
	$(F95) $(OPT) -c Mat_Mult.f90
	
BTCSmatrix.o: BTCSmatrix.f90 header.o
	$(F95) $(OPT) -c BTCSmatrix.f90
	
BTCSmatrix_PRE.o: BTCSmatrix_PRE.f90 header.o
	$(F95) $(OPT) -c BTCSmatrix_PRE.f90

norm.o: norm.f90 header.o
	$(F95) $(OPT) -c norm.f90
		
error.o: error.f90 header.o
	$(F95) $(OPT) -c error.f90	
	
error_jac.o: error_jac.f90 header.o
	$(F95) $(OPT) -c error_jac.f90	
	
jacobi.o: jacobi.f90 header.o
	$(F95) $(OPT) -c jacobi.f90		

cg.o: cg.f90 header.o
	$(F95) $(OPT) -c cg.f90

pcg.o: pcg.f90 header.o
	$(F95) $(OPT) -c pcg.f90
	
greetings.o: greetings.f90
	$(F95) $(OPT) -c greetings.f90

vecdot.o: vecdot.f90 header.o
	$(F95) $(OPT) -c vecdot.f90
	
FDsolution.o: FDsolution.f90 header.o
	$(F95) $(OPT) -c FDsolution.f90		
	
alpha.o: alpha.f90 header.o
	$(F95) $(OPT) -c alpha.f90	
	
test2D.o: test2D.f90 header.o
	$(F95) $(OPT) -c test2D.f90

test2D: test2D.o initial2D.o sparsegather.o Mat_Mult.o greetings.o BTCSmatrix.o BTCSmatrix_PRE.o norm.o error.o error_jac.o jacobi.o cg.o pcg.o vecdot.o FDsolution.o alpha.o
	 $(F95) -o test2D test2D.o header.o sparsegather.o initial2D.o greetings.o Mat_Mult.o BTCSmatrix.o BTCSmatrix_PRE.o norm.o error.o jacobi.o cg.o pcg.o vecdot.o  error_jac.o FDsolution.o alpha.o  $(LIB_LIST)

# Note that the blas are only used in the linking step
clean: 
	rm -rf *~ *.o *.mod core  
