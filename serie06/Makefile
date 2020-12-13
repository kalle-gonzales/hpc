#------------------
# Options
#------------------

# Uncomment to use BLAS routines for linear algebra.
# Make sure you have the necessary libraries installed.
#USE_BLAS = 1

# Uncomment to use OpenMP for parallelization
# and time measurement.
USE_OPENMP = 1

#------------------
# Compiler flags
#------------------

CC = gcc
RM = rm -f

CFLAGS = -Wall -g -O3 -funroll-loops -march=native -mtune=native -std=c99
LIBS   = -lm
    
ifdef USE_OPENMP
  CFLAGS += -DUSE_OPENMP -fopenmp
endif

ifdef USE_BLAS
  CFLAGS += -DUSE_BLAS
  LIBS   += -lblas -lgfortran -llapack
endif

SOURCES_LIB :=	\
	stopwatch.c \
	matrix.c \
	arithmetics.c

HEADER_LIB := \
	$(SOURCES_LIB:.c=.h)

OBJECTS_LIB := \
	$(SOURCES_LIB:.c=.o)	

SOURCES_PROGRAMS :=	\
	locks.c \
	matrixsolve.c

PROGRAMS := \
	$(SOURCES_PROGRAMS:.c=)

all: $(PROGRAMS)

$(OBJECTS_LIB): %.o: %.c
	@echo Compiling \"$<\"
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGRAMS): %: %.c $(OBJECTS_LIB)
	@echo Compiling and linking \"$<\"
	$(CC) $< $(CFLAGS) $(HEADER_LIB) -o $@ $(LDFLAGS) $(OBJECTS_LIB) $(LIBS)

clean:
	$(RM) $(PROGRAMS) $(OBJECTS_LIB)
