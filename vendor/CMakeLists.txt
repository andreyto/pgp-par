include(ExternalProject)

set(MAKEFLOW_BUILD_MPI_BACKEND FALSE CACHE BOOL "Build MPI backend for Makeflow")
if(MAKEFLOW_BUILD_MPI_BACKEND)
    set(MAKEFLOW_MPI_PATH "" CACHE PATH "MPI Home dir")
    if(NOT MAKEFLOW_MPI_PATH)
        #You can set MPI_COMPILER externally if you do not
        #want to rely on what find_package() finds.
        #find_package(MPI) will use that compiler value
        #to find derived data; otherwise it will try
        #to find the compiler first
        find_package(MPI)
        if(NOT MPI_COMPILER)
            MESSAGE(FATAL_ERROR "MPI backend for Makeflow 
            is requested by MPI, but it has not been found")
        endif()
        #mpi/bin/mpicc -> mpi/bin
        get_filename_component(mpi_path_derived "${MPI_COMPILER}" PATH)
        #mpi/bin -> mpi
        get_filename_component(mpi_path_derived "${mpi_file_derived}" PATH)
        set(MAKEFLOW_MPI_PATH ${mpi_path_derived} CACHE PATH "MPI Home dir" FORCE)
    endif()
    set(makeflow_conf_args_mpi "--with-mpi-path ${MAKEFLOW_MPI_PATH}")
endif()

set(MAKEFLOW_GLOBUS_PATH "" CACHE PATH "GLOBUS Home dir for building Makeflow")
if(MAKEFLOW_GLOBUS_PATH)
    set(makeflow_conf_args_globus "--with-globus-path ${MAKEFLOW_GLOBUS_PATH}")
endif()

ExternalProject_Add(CCTOOLS
	URL ${VENDOR_DISTRO}/cctools-3.4.3-src.at.tar.gz
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    INSTALL_DIR ${VENDOR_PREFIX_INSTALL}
    CONFIGURE_COMMAND ./configure --prefix <INSTALL_DIR>
		--without-system-sand
		--without-system-allpairs
		--without-system-wavefront
		--without-system-ftp-lite
		--without-system-s3tools
		--without-system-chirp
		--without-system-parrot
        ${makeflow_conf_args_mpi} 
        ${makeflow_conf_args_globus}
	BUILD_COMMAND make
	INSTALL_COMMAND make install
    )

#extract property from the external project and assign to
#global variable
ExternalProject_Get_Property(CCTOOLS INSTALL_DIR)
set(CCTOOLS_INSTALL_DIR "${INSTALL_DIR}" CACHE INTERNAL "")
set(MAKEFLOW "${CCTOOLS_INSTALL_DIR}/bin/makeflow" CACHE INTERNAL "")

set(pepnovo_CMakeLists "${CMAKE_CURRENT_SOURCE_DIR}/pepnovo.CMakeLists.txt")
ExternalProject_Add(PEPNOVO
    URL ${VENDOR_DISTRO}/productionPepNovo.tar.gz
    INSTALL_DIR ${VENDOR_PREFIX_INSTALL}/productionPepNovo
    #CMAKE_GENERATOR "Unix Makefiles"
    PATCH_COMMAND ${CMAKE_COMMAND} 
        -E copy 
            "${pepnovo_CMakeLists}" 
            CMakeLists.txt
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
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

