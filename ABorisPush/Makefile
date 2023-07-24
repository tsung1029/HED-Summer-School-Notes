#Makefile for 2-1/2D Electromagnetic PIC codes

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

#LEGACY =

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

#OPTS90 = -O3
OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall -std=c99

#Python
PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix
PYOPTS90 = $(OPTS90) -fPIC
PYCCOPTS = $(CCOPTS) -fPIC

LEGACY = -std=legacy

# Makefile PGI compiler with Linux

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -Wall -std=c99

#LEGACY = -Mnomain

# Makefile Nag compiler with Linux

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

#CCOPTS = -O3 -Wall -std=c99

#LEGACY =

# Makefile Nag compiler with Linux

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

#CCOPTS = -O3 -Wall -std=c99

#LEGACY =

# Makefile Cray compiler with Linux

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#LEGACY =

#

# Linkage rules

all : bpic2

python : bpic2.py

bpic2 : bpic2.o bpush2.o abpush2.o dtimer.o
	$(FC90) $(OPTS90) -o bpic2 bpic2.o wbpush2.o bpush2.o bpush2_h.o \
        wabpush2.o abpush2.o abpush2_h.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

bpush2.o : bpush2.f
	$(FC90) $(OPTS90) -o bpush2.o -c bpush2.f

bpush2_h.o : bpush2_h.f90
	$(FC90) $(OPTS90) -o bpush2_h.o -c bpush2_h.f90

abpush2.o : abpush2.f
	$(FC90) $(OPTS90) -o abpush2.o -c abpush2.f

abpush2_h.o : abpush2_h.f90
	$(FC90) $(OPTS90) -o abpush2_h.o -c abpush2_h.f90

wbpush2.o : bpush2_h.o
	$(FC90) $(OPTS90) -o wbpush2.o -c wbpush2.f90

wabpush2.o : abpush2_h.o
	$(FC90) $(OPTS90) -o wabpush2.o -c wabpush2.f90

bpic2.o : bpic2.f90 wbpush2.o wabpush2.o
	$(FC90) $(OPTS90) -o bpic2.o -c bpic2.f90

# Python

bpic2.py : bpic2lib.so dtimer.so

pybpush2_h.o : bpush2_h.f90
	$(FC90) $(PYOPTS90) -o pybpush2_h.o -c bpush2_h.f90

pyabpush2_h.o : abpush2_h.f90
	$(FC90) $(PYOPTS90) -o pyabpush2_h.o -c abpush2_h.f90

pybpush2.o : bpush2.f
	$(FC90) $(PYOPTS90) -o pybpush2.o -c bpush2.f

pyabpush2.o : abpush2.f
	$(FC90) $(PYOPTS90) -o pyabpush2.o -c abpush2.f

#double precision
bpic2lib.so : dwbpush2.f90 pybpush2_h.o pybpush2.o dwabpush2.f90 \
              pyabpush2_h.o pyabpush2.o
	$(PFC) --f90flags="$(OPTS90)" -m bpic2lib \
        -c dwbpush2.f90 dwabpush2.f90 pybpush2.o pybpush2_h.o \
        pyabpush2.o pyabpush2_h.o

#single precision
#bpic2lib.so : wbpush2.f90 pybpush2_h.o pybpush2.o wabpush2.f90 \
#              pyabpush2_h.o pyabpush2.o
#	$(PFC) --f90flags="$(OPTS90)" -m bpic2lib \
#        -c wbpush2.f90 wabpush2.f90 pybpush2.o pybpush2_h.o pyabpush2.o \
#        pyabpush2_h.o

dtimer.so : dtimer_hpy.f90 dtimer.c
	$(PFC) -m dtimer -c dtimer_hpy.f90 dtimer.c

clean :
	rm -f *.o *.mod *.pyf

clobber: clean
	rm -f bpic2 *.so

list:
	echo bpic2 python
	
