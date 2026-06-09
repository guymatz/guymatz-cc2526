#FC=g95
STDF= -std=f95
FC=gfortran
FFLAGS =

TARGET= ../bin/verlet3 
OBJECTS= verlet3.o \
         jpca15.o \
         kinds.o

%.o: %.mod

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJECTS)

verlet3.o: verlet3.f95  jpca15.o kinds.o
jpca15.o: jpca15.f kinds.o
kinds.o: kinds.f95

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95

.SUFFIXES: .f.o
.f.o:
	$(FC) $(FFLAGS) -c -I. -o $*.o $*.f

clean:
	$(RM) $(TARGET) *.mod *.o
