LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

#all: test_verlet broken_tests error_handling disabled_test
all: test_verlet test_jpca15 run

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

#STDLIB_CFLAGS := `pkg-config --cflags fortran_stdlib`
STDLIB_CFLAGS := $(shell pkg-config --cflags fortran_stdlib)
STDLIB_LIBS := $(shell pkg-config --libs fortran_stdlib)

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../src $(STDLIB_LIBS) $(STDLIB_CFLAGS)
#FFLAGS += -I../src -Wno-underflow

#test_verlet_TESTS := test_verlet.pf test_verlet_fixture.pf
test_verlet_TESTS := test_verlet.pf
test_verlet_OTHER_LIBRARIES := -L../src -lsut $(STDLIB_LIBS)
$(eval $(call make_pfunit_test,test_verlet))

test_jpca15_TESTS := test_jpca15.pf
test_jpca15_OTHER_LIBRARIES := -L../src -lsut $(STDLIB_LIBS)
$(eval $(call make_pfunit_test,test_jpca15))

#broken_tests_TESTS := test_failing.pf
#broken_tests_OTHER_LIBRARIES := -L../src -lsut
#$(eval $(call make_pfunit_test,broken_tests))
#
#error_handling_TESTS := test_error_handling.pf
#error_handling_OTHER_LIBRARIES := -L../src -lsut
#error_handling_OTHER_SRCS := throw_with_pfunit.F90
#error_handling_EXTRA_USE := throw_with_pfunit_mod
#error_handling_EXTRA_INITIALIZE := initialize_throw
#$(eval $(call make_pfunit_test,error_handling))
#
#error_handling_driver.o: throw_with_pfunit.o
#
#disabled_test_TESTS := test_disable.pf
#disabled_test_OTHER_LIBRARIES := -L../src -lsut
#$(eval $(call make_pfunit_test,disabled_test))

run:
	$(info *********   TEST OUTPUT   **************)
	./test_verlet -v
	./test_jpca15 -v


clean:
	$(RM) *.o *.mod *.a  *.inc
	$(RM) test_verlet.F90 test_failing.F90 test_error_handling.F90 test_disable.F90
	$(RM) test_verlet_fixture.F90
	$(RM) test_verlet broken_tests error_handling disabled_test
	$(RM) test_jpca15 test_jpca15.F90 
