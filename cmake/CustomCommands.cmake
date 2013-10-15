## Create cutsom_command that also sets environment variables
function(make_script_command SCRIPT_NAME SCRIPT_BODY COMMAND SCRIPT_FILE)
    set(scr_name ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT_NAME}.cmake)
    file(WRITE ${scr_name} "${SCRIPT_BODY}")
    set(cmd ${CMAKE_COMMAND} -P ${scr_name})
    set (${COMMAND} ${cmd} PARENT_SCOPE)
endfunction()

