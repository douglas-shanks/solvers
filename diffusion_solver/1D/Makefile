LIB_LIST = -llapack -lblas

# Above is the location of the blas library which is called by lublas2
# The filename with full path for this library is /usr/lib/libblas.a
# Since it is in /usr/lib we do not need to specify a path.
# The compiler looks there automatically.
F95 = gfortran
OPT = 

all:test1D

exactsoln.o: exactsoln.f90
	$(F95) -c exactsoln.f90
	
FTCS.o: FTCS.f90
	$(F95) -c FTCS.f90

FTCSmatrix.o: FTCSmatrix.f90
	$(F95) -c FTCSmatrix.f90

BTCSmatrix.o: BTCSmatrix.f90
	$(F95) -c BTCSmatrix.f90

BTCS_rhs.o: BTCS_rhs.f90
	$(F95) -c BTCS_rhs.f90

CNmatrix.o: CNmatrix.f90
	$(F95) -c CNmatrix.f90

CN_rhs.o: CN_rhs.f90
	$(F95) -c CN_rhs.f90
	
test1D.o: test1D.f90
	$(F95) -c test1D.f90

norm.o: norm.f90
	$(F95) -c norm.f90
	
test1D: test1D.o exactsoln.o FTCS.o FTCSmatrix.o BTCSmatrix.o CNmatrix.o BTCS_rhs.o CN_rhs.o norm.o
	 $(F95) -o test1D test1D.o exactsoln.o FTCS.o FTCSmatrix.o BTCSmatrix.o CNmatrix.o BTCS_rhs.o CN_rhs.o norm.o  $(LIB_LIST)

# Note that the blas are only used in the linking step
clean: 
	rm -rf *~ *.o *.mod core  
