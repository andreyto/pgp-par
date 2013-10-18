include(Checks)
include(CMakeParseArguments)

## Write CMake script_body to a file and return command line
## for running the script and the file name of the script
## These values can be used as arguments of custom_command():
## add_custom_command(COMMAND script_command DEPENDS script_file ...).
## If your script body is a multiline string that contains:
## "
## set(ENV{VAR}=${val})
## execute_process(COMMAND whatever ...)
## "
## and 'val' is defined in your CMakeLists.txt,
## then this is the way to generate custom_command target that 
## is built under a modified shell environment.
function(make_script_command) 
    set(oneValueArgs SCRIPT_NAME SCRIPT_BODY VAR_SCRIPT_COMMAND VAR_SCRIPT_FILE)
    cmake_parse_arguments(arg "" "${oneValueArgs}" "" ${ARGN})
    at_assert_defined_prefixed(arg ${oneValueArgs})
    set(scr_name ${CMAKE_CURRENT_BINARY_DIR}/${arg_SCRIPT_NAME}.cmake)
    file(WRITE ${scr_name} "${arg_SCRIPT_BODY}")
    set(cmd ${CMAKE_COMMAND} -P ${scr_name})
    set (${arg_VAR_SCRIPT_COMMAND} ${cmd} PARENT_SCOPE)
    set (${arg_VAR_SCRIPT_FILE} ${scr_name} PARENT_SCOPE)
endfunction()

