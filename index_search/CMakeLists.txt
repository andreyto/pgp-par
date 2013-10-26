## Python extension and dependencies (Biopython)
include(CustomCommands)

set(PY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(PY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/py_build")
set(SETUP_PY_IN "${PY_SOURCE_DIR}/setup.py.in")
#easy_install needs the source dir as the last argument,
#and expects that dir to contain setup.py,
#plus it will download distribute tarball into it.
#We can either copy all python source into cmake binary dir
#and build in that copy, or build in the original source dir.
#We choose the late so that the development cycle modify/test 
#is easier.
set(SETUP_PY    "${PY_SOURCE_DIR}/setup.py")
set(DEPS        "${PY_SOURCE_DIR}/IndexSearch/__init__.py")
set(OUTPUT      "${PY_BUILD_DIR}/timestamp")

set(PY_BUILD_INCLUDE_DIRS ${SEQAN_INCLUDE_DIRS})

configure_file(${SETUP_PY_IN} ${SETUP_PY} @ONLY)

file(MAKE_DIRECTORY ${PY_BUILD_INSTALL_PREFIX})
file(MAKE_DIRECTORY ${PY_BUILD_DIR})

#Doing just 'setup.py build' or 'setup.py develop'
#does not generate
#entry point scripts, thus making in-build testing
#impossible. Thus, we do 'setup.py install' twice -
#during build and during install, into different locations.

#We use --site-dirs argument for easy_install to force installation
#into a directory which is not in the PYTHONPATH 
#(i.e. .pth files would not work where). easy_install still complains,
#in a slightly differnt way, so we also modify the setup.py above.

#If we ever need to modify the enivronment for custom_command,
#that can be only done by running a generated cmake script
#as the custom command, such as:
#add_custom_command( 
#  COMMAND ${CMAKE_COMMAND} -DVAR=VAL 
#    -P ${PROJECT_BINARY_DIR}/script.cmake 
#  DEPENDS ${PROJECT_BINARY_DIR}/script.cmake 
#  WORKING_DIRECTORY ${PROJECT_BINARY_DIR} 
#  COMMENT "Running script.cmake" 
#  VERBATIM) 
#configure_file(myscript.cmake.in myscript.cmake @ONLY)
#myscript.cmake:
#set(ENV{PATH} ...)
#execute_process(....)
#
#The reason for the above is because custom_command will
#be executed at build time (e.g. when 'make' is running),
#so any changes we make by set(ENV{}) here will be already lost.
#
#For install, use install(SCRIPT ...) signature.


set(cmd_build_install "${PYTHON} ${SETUP_PY} easy_install 
					--install-dir ${PY_BUILD_INSTALL_PREFIX}
					--site-dirs ${PY_BUILD_INSTALL_PREFIX}
					--build-directory ${PY_BUILD_DIR} 
					--always-copy 
					.")


set(cmd_develop ${PYTHON} ${SETUP_PY} develop 
					--install-dir ${PY_BUILD_INSTALL_PREFIX} 
					--build-directory ${PY_BUILD_DIR})

set(sitecustomize_body "
##This file causes .pth files to be processed if the directory
##of this file is on PYTHONPATH. Placing this file in such a directory
##is the way to create an additional site-packages like directory in which
##eggs can be installed without having to edit any system or user wide Python
##config files.
import site
import os,inspect
site.addsitedir(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
")

set(cmd "
set(ENV{CFLAGS} -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0)
set(ENV{CXXFLAGS} -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0)
execute_process(
    COMMAND ${cmd_build_install}
    RESULT_VARIABLE cmd_error_status
    )
if( cmd_error_status )
    message(FATAL_ERROR 'execute_process failed')
endif()
file(WRITE ${PY_BUILD_INSTALL_PREFIX}/sitecustomize.py \"${sitecustomize_body}\")
")

make_script_command(SCRIPT_NAME py_build 
    SCRIPT_BODY "${cmd}" 
    VAR_SCRIPT_COMMAND cmd_build_install_scr 
    VAR_SCRIPT_FILE build_install_scr)

add_custom_command(OUTPUT ${OUTPUT}
                    COMMAND ${cmd_build_install_scr}
                    #COMMAND ${cmd_develop}
                    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                    DEPENDS ${build_install_scr} ${DEPS}
					WORKING_DIRECTORY ${PY_SOURCE_DIR} 
					VERBATIM
					COMMENT "Command to build and install Python modules into build dir")

add_custom_target(PY_BUILD ALL 
	DEPENDS ${OUTPUT} 
	COMMENT "Target to build and install Python modules into build dir")

set(cmd_install ${PYTHON} ${SETUP_PY} easy_install 
					--install-dir ${PY_INSTALL_PREFIX} 
					--build-directory ${PY_BUILD_DIR} 
					--site-dirs ${PY_INSTALL_PREFIX} 
					--always-copy .)

#install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PY_INSTALL_PREFIX} ERROR_QUIET ERROR_VARIABLE err)")
#install(CODE "execute_process(COMMAND ${cmd_install} WORKING_DIRECTORY ${PY_SOURCE_DIR})")
#trailing slash copies files but not the directory itself
install(DIRECTORY ${PY_BUILD_INSTALL_PREFIX}/ DESTINATION ${PY_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)

