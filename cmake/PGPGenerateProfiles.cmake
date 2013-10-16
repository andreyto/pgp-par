
set(config_dir "${CMAKE_INSTALL_PREFIX}/config")
set(bin_dir "${CMAKE_INSTALL_PREFIX}/bin")

set(PGP_LOGIN_RC "${config_dir}/pgp_login.rc")
set(PGP_TASK_RC "${config_dir}/pgp_task.rc")
set(PGP_PREPARE "${bin_dir}/pgp_prepare")
set(PGP_PREPARE_SUBMIT "${bin_dir}/pgp_prepare_submit")
set(PGP_RUN_SUBMIT "${bin_dir}/pgp_run_submit")
set(PGP_PREPARE_QSUB_IN "${config_dir}/pgp_prepare.qsub")
set(PGP_RUN_QSUB_IN "${config_dir}/pgp_run.qsub")
set(PGP_CONFIG "${config_dir}/pgp.ini")
set(PGP_ROOT "${CMAKE_INSTALL_PREFIX}")
set(PGP_HOME "${PGP_INSTALL_PREFIX}")

foreach(dir_name in config bin)
    set(dir_out "${CMAKE_BINARY_DIR}/${dir_name}")
    file(MAKE_DIRECTORY "${dir_out}")
    file(GLOB arch_files_in ${PROJECT_SOURCE_DIR}/${dir_name}/${PGP_TARGET_ENV}/*)
    file(GLOB noarch_files_in ${PROJECT_SOURCE_DIR}/${dir_name}/noarch/*)
    foreach( file_in IN LISTS arch_files_in noarch_files_in)
        get_filename_component(file_out "${file_in}" NAME)
        string(REGEX REPLACE ".in\$" "" file_out "${file_out}")
        set(file_out "${dir_out}/${file_out}")
        configure_file(${file_in} ${file_out} @ONLY)
    endforeach()
endforeach()

install(DIRECTORY "${CMAKE_BINARY_DIR}/bin" 
    DESTINATION "${CMAKE_INSTALL_PREFIX}" 
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE GROUP_READ)

install(DIRECTORY "${CMAKE_BINARY_DIR}/config" 
    DESTINATION "${CMAKE_INSTALL_PREFIX}")

