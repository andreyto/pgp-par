#set(CMAKE_CXX_FLAGS "-O3" CACHE INTERNAL "")

set(MAKEFLOW_BUILD_MPI_BACKEND TRUE CACHE BOOL "Build MPI backend for Makeflow")
set(MAKEFLOW_MPI_PATH "$ENV{MPICH_HOME}" CACHE PATH "MPI Home dir for building Makeflow")
set(MAKEFLOW_GLOBUS_PATH "$ENV{GLOBUS_PATH}" CACHE PATH "GLOBUS Home dir for building Makeflow")

