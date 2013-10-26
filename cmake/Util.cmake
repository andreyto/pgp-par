## Get delimiter for entries within PATH, PYTHONPATH etc
function(get_path_separator python var_path_sep)
    execute_process(COMMAND ${python} -c "'import os; print os.pathsep()'"
        RESULT_VARIABLE cmd_error_status
        OUTPUT_VARIABLE path_sep)
    if( cmd_error_status )
        message(FATAL_ERROR execute_process failed)
    endif()
    set(${var_path_sep} ${path_sep} PARENT_SCOPE)
endfunction()

