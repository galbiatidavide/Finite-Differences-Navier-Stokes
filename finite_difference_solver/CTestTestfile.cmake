# CMake generated Testfile for 
# Source directory: /home/dave/lifex/tests/finite_difference_solver
# Build directory: /home/dave/lifex/tests/finite_difference_solver
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_finite_differences "/u/sw/toolchains/gcc-glibc/11.2.0/base/bin/mpiexec" "-n" "4" "/home/dave/lifex/tests/finite_difference_solver/lifex_test_finite_differences")
set_tests_properties(test_finite_differences PROPERTIES  LABELS "test_hard" WORKING_DIRECTORY "/home/dave/lifex/tests/finite_difference_solver" _BACKTRACE_TRIPLES "/home/dave/lifex/cmake/setup_lifex_tests.cmake;83;add_test;/home/dave/lifex/cmake/setup_lifex_tests.cmake;113;_lifex_add_test;/home/dave/lifex/tests/finite_difference_solver/CMakeLists.txt;23;lifex_add_test_hard_mpi;/home/dave/lifex/tests/finite_difference_solver/CMakeLists.txt;0;")
