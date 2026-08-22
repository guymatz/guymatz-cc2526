# Necessary if the installation directory is not in PKG_CONFIG_PATH
#install_dir := path/to/install_dir
#export PKG_CONFIG_PATH := $(install_dir)/lib/pkgconfig:$(PKG_CONFIG_PATH)

# On some OS the pkgconfig file could be installed in lib64 instead of lib
# export PKG_CONFIG_PATH := $(install_dir)/lib64/pkgconfig:$(PKG_CONFIG_PATH)

STDLIB_CFLAGS := `pkg-config --cflags fortran_stdlib`
STDLIB_LIBS := `pkg-config --libs fortran_stdlib`


#FC=g95
# THIS BREAKS LINKING
#STDF= -std=f95
FC=gfortran
#FFLAGS = -g -fbacktrace -fcheck=all -Wall
FFLAGS = -g -fbacktrace -Wall

PFUNIT_DIR = /usr/local/pfunit

LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

#all: clean $(BIN_TARGET) tests

BIN_TARGET = ../bin/verlet3
#TEST_TARGET = ../bin/test-verlet
BIN_OBJECTS= main.o \
         verlet.o \
         jpca15.o \
         kinds.o

#TEST_OBJECTS= test-verlet.o \
#         verlet.o \
#         jpca15.o \
#         kinds.o
#         # verlet.o \

%.o: %.mod

$(BIN_TARGET): libsut.a $(BIN_OBJECTS)
	$(FC) $(FFLAGS) -o $(BIN_TARGET) $(BIN_OBJECTS) $(STDLIB_LIBS)

#$(TEST_TARGET): $(TEST_OBJECTS)
#	$(FC) $(FFLAGS) -o $(TEST_TARGET) $(TEST_OBJECTS)

main.o: main.F90  jpca15.o kinds.o verlet.o
#test-verlet.o: test-verlet.f95  jpca15.o kinds.o verlet.o
verlet.o: verlet.f95 kinds.o
jpca15.o: jpca15.f kinds.o
kinds.o: kinds.f95

all: clean $(BIN_TARGET) run

#.PHONY: clean $(BIN_TARGET) $(TEST_TARGET) run

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95 $(STDLIB_CFLAGS) $(STDLIB_LIBS)

.SUFFIXES: .F90 .o
.F90.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.F90 $(STDLIB_CFLAGS) $(STDLIB_LIBS)

.SUFFIXES: .f.o
.f.o:
	$(FC) $(FFLAGS) -c -I. -o $*.o $*.f $(STDLIB_CFLAGS) $(STDLIB_LIBS)

libsut.a: $(BIN_OBJECTS) kinds.o 
	$(AR) -r $@ $?

clean:
	$(RM) $(BIN_TARGET) $(TEST_TARGET) libsut.a *.mod *.o tests tests.inc test_verlet.F90

run: clean $(BIN_TARGET)
	$(TBIN_ARGET)

#test: clean $(TEST_TARGET)
#	$(TEST_TARGET)

tests: libsut.a

tests_TESTS := test_verlet.pf
tests_REGISTRY :=
tests_OTHER_SOURCES :=
tests_OTHER_LIBRARIES := -L. -lsut
tests_OTHER_INCS :=

$(eval $(call make_pfunit_test,tests))
