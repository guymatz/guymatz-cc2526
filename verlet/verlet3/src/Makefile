#FC=g95
# THIS BREAKS LINKING
#STDF= -std=f95
FC=gfortran
FFLAGS =

TARGET= ../bin/verlet3 
OBJECTS= main.o \
         jpca15.o \
         kinds.o
         # verlet.o \

%.o: %.mod

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJECTS)

main.o: main.f95  jpca15.o kinds.o# verlet.o
# verlet.o: verlet.f95 kinds.o
jpca15.o: jpca15.f kinds.o
kinds.o: kinds.f95

.PHONY: clean $(TARGET) run

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95

.SUFFIXES: .f.o
.f.o:
	$(FC) $(FFLAGS) -c -I. -o $*.o $*.f

clean:
	$(RM) $(TARGET) *.mod *.o

run: clean $(TARGET)
	$(TARGET)
