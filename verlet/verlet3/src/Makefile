#FC=g95
# THIS BREAKS LINKING
#STDF= -std=f95
FC=gfortran
#FFLAGS = -g -fbacktrace -fcheck=all -Wall
FFLAGS = -g -fbacktrace -Wall

PFUNIT_DIR = /usr/local/pfunit

include $(PFUNIT_DIR)/PFUNIT/include/PFUNIT.mk
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

# all: TARGET1 TARGET2

BIN_TARGET = ../bin/verlet3
TEST_TARGET = ../bin/test-verlet
BIN_OBJECTS= main.o \
         verlet.o \
         jpca15.o \
         kinds.o

TEST_OBJECTS= test-verlet.o \
         verlet.o \
         jpca15.o \
         kinds.o
         # verlet.o \

%.o: %.mod

$(BIN_TARGET): $(BIN_OBJECTS)
	$(FC) $(FFLAGS) -o $(BIN_TARGET) $(BIN_OBJECTS)

$(TEST_TARGET): $(TEST_OBJECTS)
	$(FC) $(FFLAGS) -o $(TEST_TARGET) $(TEST_OBJECTS)

main.o: main.F90  jpca15.o kinds.o verlet.o
test-verlet.o: test-verlet.f95  jpca15.o kinds.o verlet.o
verlet.o: verlet.f95 kinds.o
jpca15.o: jpca15.f kinds.o
kinds.o: kinds.f95

.PHONY: clean $(BIN_TARGET) $(TEST_TARGET) run

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95

.SUFFIXES: .F90 .o
.F90.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.F90

.SUFFIXES: .f.o
.f.o:
	$(FC) $(FFLAGS) -c -I. -o $*.o $*.f

clean:
	$(RM) $(BIN_TARGET) $(TEST_TARGET) *.mod *.o

run: clean $(BIN_TARGET)
	$(TBIN_ARGET)

test: clean $(TEST_TARGET)
	$(TEST_TARGET)
