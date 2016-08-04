LIB_LIST = -llapack -lblas

# Above is the location of the blas library which is called by lublas2
# The filename with full path for this library is /usr/lib/libblas.a
# Since it is in /usr/lib we do not need to specify a path.
# The compiler looks there automatically.
F95 = gfortran
OPT = -O

all:test2D

header.o: header.f90
	$(F95) $(OPT) -c header.f90
	
exactsoln2D.o: exactsoln2D.f90 header.o
	$(F95) $(OPT) -c exactsoln2D.f90	
	
Mat_Mult.o: Mat_Mult.f90 header.o
	$(F95) $(OPT) -c Mat_Mult.f90
	
FTCSmatrix.o: FTCSmatrix.f90 header.o
	$(F95) $(OPT) -c FTCSmatrix.f90
	
BTCSmatrix.o: BTCSmatrix.f90 header.o
	$(F95) $(OPT) -c BTCSmatrix.f90
	
BTCSmatrix_PRE.o: BTCSmatrix_PRE.f90 header.o
	$(F95) $(OPT) -c BTCSmatrix_PRE.f90
	
CNmatrix.o: CNmatrix.f90 header.o
	$(F95) $(OPT) -c CNmatrix.f90	

CNmatrix_PRE.o: CNmatrix_PRE.f90 header.o
	$(F95) $(OPT) -c CNmatrix_PRE.f90	

CNrhs.o: CNrhs.f90 header.o
	$(F95) $(OPT) -c CNrhs.f90

norm.o: norm.f90 header.o
	$(F95) $(OPT) -c norm.f90
		
error.o: error.f90 header.o
	$(F95) $(OPT) -c error.f90	
	
jacobi.o: jacobi.f90 header.o
	$(F95) $(OPT) -c jacobi.f90		

cg.o: cg.f90 header.o
	$(F95) $(OPT) -c cg.f90

pcg.o: pcg.f90 header.o
	$(F95) $(OPT) -c pcg.f90

vecdot.o: vecdot.f90 header.o
	$(F95) $(OPT) -c vecdot.f90
	
FDsolution.o: FDsolution.f90 header.o
	$(F95) $(OPT) -c FDsolution.f90		

test2D.o: test2D.f90 header.o
	$(F95) $(OPT) -c test2D.f90

test2D: test2D.o exactsoln2D.o Mat_Mult.o FTCSmatrix.o BTCSmatrix.o BTCSmatrix_PRE.o CNmatrix.o CNmatrix_PRE.o CNrhs.o norm.o error.o jacobi.o cg.o pcg.o vecdot.o FDsolution.o
	 $(F95) -o test2D test2D.o header.o exactsoln2D.o Mat_Mult.o FTCSmatrix.o BTCSmatrix.o BTCSmatrix_PRE.o CNmatrix.o CNmatrix_PRE.o CNrhs.o norm.o error.o jacobi.o cg.o pcg.o vecdot.o  FDsolution.o  $(LIB_LIST)

# Note that the blas are only used in the linking step
clean: 
	rm -rf *~ *.o *.mod core  
