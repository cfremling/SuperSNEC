# some definitions controlling compilation
export F90=gfortran
#export F90FLAGS=-O3 -mcpu=native -funroll-loops -ffast-math
export F90FLAGS=-Ofast -fprotect-parens -mcpu=native -funroll-loops -flto
export LDFLAGS=$(F90FLAGS) 
export LAPACKLIBS=-framework Accelerate
export LIBS+=$(LAPACKLIBS)
