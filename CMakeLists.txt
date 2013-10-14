cmake_minimum_required(VERSION 2.8.0)
#cmake_minimum_required(VERSION 2.8.5)

PROJECT(PGP)

SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 2 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )

SET( ${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

set(PACKAGE_VERSION ${${PROJECT_NAME}_VERSION})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#MESSAGE("CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")

# For make-based builds, defines make target named test.
# For Visual Studio builds, defines Visual Studio project named RUN_TESTS.
enable_testing()

set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/test_data
    CACHE FILEPATH "Root directory for input data files used by the test cases")

find_program(PYTHON "python")

set(VENDOR_DISTRO ${PROJECT_SOURCE_DIR}/vendor)
set(VENDOR_PREFIX_INSTALL ${CMAKE_INSTALL_PREFIX}/vendor)

#We want to be able to test from the build diretory,
#without having to run install every time. We unify
#the build output locations below.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(PY_INSTALL_SUBDIR py)
set(PY_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_SUBDIR}")
set(PY_BUILD_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${PY_INSTALL_SUBDIR}")

set(PGP_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/proteogenomics")

set(WRAPPER_SCRIPT_ROOT pgp)

#set(WRAPPER_SCRIPT_IN "${CMAKE_SOURCE_DIR}/bin/${WRAPPER_SCRIPT_ROOT}.in")
#set(WRAPPER_SCRIPT "${CMAKE_BINARY_DIR}/${WRAPPER_SCRIPT_ROOT}")

#configure_file(${WRAPPER_SCRIPT_IN} ${WRAPPER_SCRIPT} @ONLY)
#install(PROGRAMS ${WRAPPER_SCRIPT} DESTINATION ${CMAKE_INSTALL_PREFIX})

add_subdirectory(vendor)
add_subdirectory(proteogenomics)
add_subdirectory(index_search)


