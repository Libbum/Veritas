cmake_minimum_required (VERSION 3.1 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 14)
project (Veritas)
set (Veritas_VERSION_MAJOR 1)
set (Veritas_VERSION_MINOR 0)
set (Veritas_VERSION_REVISION 0)

set (CMAKE_BUILD_TYPE Release)
option (USE_MKL "Use MKL Library" ON)

set (CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/Modules/" CACHE PATH "CMake module path" FORCE)

find_package (OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if (USE_MKL)
   find_package (MKL)
   if (MKL_FOUND)
      #If we use the intel compilers, the CMAKE_CXX_STANDARD variable isn't read. Also we add a few more optimisations.
      if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
         #Default IPO, drop superflous warnings
         set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -parallel -ansi-alias -DMKL_ILP64 -I${MKL_ROOT}/include -std=c++14 -fp-model fast -wd11074 -wd11076")
      else()
         #Assuming GNU
         set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64 -m64 -I${MKL_ROOT}/include")
      endif()
      add_definitions (-DUSINGMKL)
   else()
      message ("-- MKL Library Not Found. Falling back to LAPACK.")
      find_package (LAPACK REQUIRED)
   endif()
else()
   find_package (LAPACK REQUIRED)
endif()



include (GetGitRevisionDescription)
get_git_head_revision (GIT_REFSPEC GIT_SHA1)
configure_file (gitsha1.hpp.in ${CMAKE_BINARY_DIR}/generated/gitsha1.hpp @ONLY)
include_directories (${CMAKE_BINARY_DIR}/generated)

add_executable (veritas BoundaryCondition.cpp Level.cpp Rectangle.cpp SolverManager.cpp EMSolver.cpp Mesh.cpp Settings.cpp veritas.cpp)
if (MKL_FOUND)
   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      target_link_libraries (veritas -Wl,--start-group ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKL_ROOT}/lib/intel64/libmkl_core.a ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl)
   else()
      target_link_libraries (veritas -Wl,--start-group ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKL_ROOT}/lib/intel64/libmkl_core.a ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl)
   endif()
else()
   target_link_libraries (veritas lapack blas)
endif()
set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
install (TARGETS veritas DESTINATION bin)
install (DIRECTORY DESTINATION bin/output/rectangleData)
