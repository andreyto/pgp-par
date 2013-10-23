set(test_archive "${TEST_DATA_DIR}/Bacillus.anthracis.sterne.PNNL.chunk.4sp.tar.gz")
add_test(NAME test_small 
    COMMAND ${PGP_WRAPPER} ${PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/test_run.py
        ${test_archive} 
        ${CMAKE_INSTALL_PREFIX}
        ${PGP_TEST_PGP_RUN_EXTRA_OPTIONS})

