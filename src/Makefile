
# Set MKLROOT if it is not defined
#MKLROOT = /opt/intel/oneapi/mkl/latest

# Compilers
CC = gcc
FC = gfortran

OBJDIR = ../obj/
BINDIR = ../bin/

# C object files
C_OBJS = $(OBJDIR)bound.o $(OBJDIR)bounding_procedure.o $(OBJDIR)build_subproblem_data.o \
         	$(OBJDIR)miscellaneous_functions.o $(OBJDIR)instance_data.o \
         	$(OBJDIR)b_and_b_framework.o $(OBJDIR)srflp.o $(OBJDIR)heuristics.o \
         	$(OBJDIR)initial_heuristics.o $(OBJDIR)inequalities.o \
         	$(OBJDIR)parameters.o $(OBJDIR)branching.o

# Fortran object files
FORTRAN_OBJS = $(OBJDIR)blas.o $(OBJDIR)lbfgsb.o $(OBJDIR)timer.o $(OBJDIR)linpack.o $(OBJDIR)redirect_output.o

# Determine the machine architecture
ifeq ($(shell uname -m), x86_64)
	ARCH = -m64
else
	ARCH = -m32
endif

MODE = -O3 -DNDEBUG
#MODE = -g3

LINALG = -m64  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

CFLAGS = -Wall $(MODE) $(ARCH) -lpthread -I${MKLROOT}/include

FFLAGS = -frecursive -fomit-frame-pointer -fPIC $(MODE) $(ARCH)

### Rules ###

.PHONY : all clean

# Default rule
all : $(BINDIR)srflp_bb

# Rules for the final binary
$(BINDIR)srflp_bb : $(C_OBJS) $(FORTRAN_OBJS)
	$(FC) -o $@ $^ $(CFLAGS) $(LINALG)

# Compile all C files
$(OBJDIR)%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@
	
# Compile all Fortran files
$(OBJDIR)%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

clean :
	rm $(BINDIR)srflp_bb $(C_OBJS) $(FORTRAN_OBJS)
