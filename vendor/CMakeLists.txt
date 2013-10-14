include(ExternalProject)

ExternalProject_Add(CCTOOLS
	URL ${VENDOR_DISTRO}/cctools-3.4.3-src.at.tar.gz
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    INSTALL_DIR ${VENDOR_PREFIX_INSTALL}
    CONFIGURE_COMMAND ./configure --prefix <INSTALL_DIR>
	BUILD_COMMAND make
	INSTALL_COMMAND make install
    )

#extract property from the external project and assign to
#global variable
ExternalProject_Get_Property(CCTOOLS INSTALL_DIR)
set(CCTOOLS_INSTALL_DIR ${INSTALL_DIR} CACHE INTERNAL "")

#TODO: pepnovo uses g++ hardcoded; convert its build to CMake
ExternalProject_Add(PEPNOVO
    URL ${VENDOR_DISTRO}/productionPepNovo.tar.gz
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    INSTALL_DIR ${VENDOR_PREFIX_INSTALL}/productionPepNovo
    CONFIGURE_COMMAND ""
	BUILD_COMMAND cd src && make
    INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR> && 
    cd <SOURCE_DIR> && cp -a copyright.txt  Models  
    MS-Filter.readme.txt readme.txt <INSTALL_DIR> && 
    cp -a src/PepNovo_bin <INSTALL_DIR>
    )

ExternalProject_Get_Property(PEPNOVO INSTALL_DIR)
set(PEPNOVO_INSTALL_DIR ${INSTALL_DIR} CACHE INTERNAL "")

## It's a jar, no build
ExternalProject_Add(MSGF
    URL ${VENDOR_DISTRO}/MSGF.tar.gz
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    INSTALL_DIR ${VENDOR_PREFIX_INSTALL}/MSGF
    CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR> &&
    ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR> <INSTALL_DIR>  
    #cp -a  <SOURCE_DIR>/* <INSTALL_DIR>
    )

ExternalProject_Get_Property(MSGF INSTALL_DIR)
set(MSGF_INSTALL_DIR ${INSTALL_DIR} CACHE INTERNAL "")

## We use it as header only lib during build,
## no need to install
ExternalProject_Add(SEQAN
    URL ${VENDOR_DISTRO}/seqan-1.3.1.tar.gz
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
    INSTALL_COMMAND ""
    )

ExternalProject_Get_Property(SEQAN SOURCE_DIR)
set(SEQAN_INCLUDE_DIRS ${SOURCE_DIR} CACHE INTERNAL "")

set(VENDOR_ARCHIVE_BIOPYTHON ${VENDOR_DISTRO}/biopython-1.62.at.tar.gz CACHE INTERNAL "")

