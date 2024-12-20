# CMake generated Testfile for 
# Source directory: /home/dave/new_lifex/lifex/tests/finite_difference_solver
# Build directory: /home/dave/new_lifex/lifex/tests/finite_difference_solver
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_stokes_finale "/home/dave/Downloads/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/mpiexec" "-n" "4" "/home/dave/new_lifex/lifex/tests/finite_difference_solver/lifex_test_stokes_finale")
set_tests_properties(test_stokes_finale PROPERTIES  LABELS "test_hard" WORKING_DIRECTORY "/home/dave/new_lifex/lifex/tests/finite_difference_solver" _BACKTRACE_TRIPLES "/home/dave/new_lifex/lifex/cmake/setup_lifex_tests.cmake;83;add_test;/home/dave/new_lifex/lifex/cmake/setup_lifex_tests.cmake;113;_lifex_add_test;/home/dave/new_lifex/lifex/tests/finite_difference_solver/CMakeLists.txt;23;lifex_add_test_hard_mpi;/home/dave/new_lifex/lifex/tests/finite_difference_solver/CMakeLists.txt;0;")
