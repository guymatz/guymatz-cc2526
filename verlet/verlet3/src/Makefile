#FC=g95
# THIS BREAKS LINKING
#STDF= -std=f95
FC=gfortran
#FFLAGS = -g -fbacktrace -fcheck=all -Wall
FFLAGS = -g -fbacktrace -Wall

# all: TARGET1 TARGET2

TARGET1 = ../bin/verlet3
TARGET2 = ../bin/test-verlet
OBJECTS1= main.o \
         verlet.o \
         jpca15.o \
         kinds.o

OBJECTS2= test-verlet.o \
         verlet.o \
         jpca15.o \
         kinds.o
         # verlet.o \

%.o: %.mod

$(TARGET1): $(OBJECTS1)
	$(FC) $(FFLAGS) -o $(TARGET1) $(OBJECTS1)

$(TARGET2): $(OBJECTS2)
	$(FC) $(FFLAGS) -o $(TARGET2) $(OBJECTS2)

main.o: main.f95  jpca15.o kinds.o verlet.o
test-verlet.o: test-verlet.f95  jpca15.o kinds.o verlet.o
verlet.o: verlet.f95 kinds.o
jpca15.o: jpca15.f kinds.o
kinds.o: kinds.f95

.PHONY: clean $(TARGET1) $(TARGET2) run

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95

.SUFFIXES: .f.o
.f.o:
	$(FC) $(FFLAGS) -c -I. -o $*.o $*.f

clean:
	$(RM) $(TARGET1) $(TARGET2) *.mod *.o

run: clean $(TARGET1)
	$(TARGET1)

test: clean $(TARGET2)
	$(TARGET2)
