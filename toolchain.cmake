#set(CMAKE_CXX_FLAGS "/DBOOST_ALL_DYN_LINK /D_SCL_SECURE_NO_WARNINGS /bigobj /EHsc /D_HDF5USEDLL_ /DWIN32 /DGSL_DLL" CACHE INTERNAL "")

#set(CMAKE_PREFIX_PATH
#	${PRODDL_DEPS_HOME}
#	${FFTW3_DIR} CACHE INTERNAL "")
#set(ENV{PATH} "$ENV{PATH};${PRODDL_DEPS_HOME}/bin;${SWIG_DIR};${HDF5_DIR}/bin" CACHE INTERNAL "")
#set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" CACHE INTERNAL "")
#include_directories(${PRODDL_DEPS_HOME}/include ${PRODDL_DEPS_HOME}/include/opencv ${PRODDL_DEPS_HOME}/include/opencv2)
#link_directories(${PRODDL_DEPS_HOME}/bin)
#include(${OpenCV_DIR}/OpenCVConfig.cmake)

set(CMAKE_INSTALL_PREFIX /home/atovtchi/work/PGP.install CACHE INTERNAL "")
set(CMAKE_VERBOSE_MAKEFILE ON CACHE INTERNAL "")