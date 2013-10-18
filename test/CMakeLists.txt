
set(test_archive "${TEST_DATA_DIR}/Bacillus.anthracis.sterne.PNNL.chunk.4sp.tar.gz")
add_test(NAME test_small 
    COMMAND ${PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/test_run.py 
        --test-archive ${test_archive} 
        --install-prefix ${CMAKE_INSTALL_PREFIX}
        --test-source-dir ${CMAKE_CURRENT_SOURCE_DIR})

