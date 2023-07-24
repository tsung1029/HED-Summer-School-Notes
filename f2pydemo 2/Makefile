#-----------------------------------------------------------------------
#Makefile for f2pydemo codes: f2pydemo1,f2pydemo2,plasmalib,rplasmalib

FC90 = gfortran
CC = gcc

PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix
PYOPTS90 = $(OPTS90) -fPIC
PYCCOPTS = $(CCOPTS) -fPIC

all: f2pydemo1 f2pydemo2

python : f2pydemo1.py f2pydemo2.py

f2pydemo1 : f2pydemo1.o pushlib1.o
	$(FC90) $(OPTS90) -o f2pydemo1 f2pydemo1.o initlib1.o pushmod1.o \
       pushlib1.o

f2pydemo2 : f2pydemo2.o initlib1.o rpushlib1.o
	$(FC90) $(OPTS90) -o f2pydemo2 f2pydemo2.o initlib1.o rpushlib1.o

pushlib1.o : pushlib1.f
	$(FC90) $(OPTS90) -o pushlib1.o -c pushlib1.f

initlib1.o : initlib1.f90
	$(FC90) $(OPTS90) -o initlib1.o -c initlib1.f90

rpushlib1.o : rpushlib1.f90
	$(FC90) $(OPTS90) -o rpushlib1.o -c rpushlib1.f90

pushmod1.o : pushmod1.f90
	$(FC90) $(OPTS90) -o pushmod1.o -c pushmod1.f90

rpushmod1.o : rpushmod1.f90 rpushlib1.o
	$(FC90) $(OPTS90) -o rpushmod1.o -c rpushmod1.f90

f2pydemo1.o : f2pydemo1.f90 initlib1.o pushmod1.o
	$(FC90) $(OPTS90) -o f2pydemo1.o -c f2pydemo1.f90

f2pydemo2.o : f2pydemo2.f90 initlib1.o rpushlib1.o
	$(FC90) $(OPTS90) -o f2pydemo2.o -c f2pydemo2.f90

# Python

f2pydemo1.py : plasmalib.so

f2pydemo2.py : rplasmalib.so

pypushlib1.o : pushlib1.f
	$(FC90) $(PYOPTS90) -o pypushlib1.o -c pushlib1.f

pyrpushlib1.o : rpushlib1.f90
	$(FC90) $(PYOPTS90) -o pyrpushlib1.o -c rpushlib1.f90

plasmalib.so : pypushlib1.o
	$(PFC) --f90flags="$(OPTS90)" -m plasmalib -c initlib1.f90 \
	   pushmod1.f90 pypushlib1.o

rplasmalib.so : pyrpushlib1.o
	$(PFC) --f90flags="$(OPTS90)" -m rplasmalib -c initlib1.f90 \
	   rpushmod1.f90 pyrpushlib1.o

clean: 
	rm -f *~ *.mod *.o

clobber: clean
	rm -f f2pydemo1 f2pydemo2 *.so
