#FC=g95
#STDF= -std=F
FC=gfortran
FFLAGS =

TARGET= ../bin/verlet3 
OBJECTS= verlet3.o \
         kinds.o \
         jpca15.o

%.o: %.mod

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJECTS)

verlet3.o: verlet3.f95 kinds.o jpca15.o
jpca15.o: jpca15.f95 kinds.o
kinds.o: kinds.f95

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95

clean:
	$(RM) $(TARGET) *.mod *.o
