##CMake is quite difficult (if not impossible) when you need to pass a shell-quoted string
##like below to some subprocess unchanged. The quotations and STRING variable
##type below do work (with extra shlex splitting done in the test Python code).
##Follow this example when editing the value, or you will break the test.
set(PGP_TEST_PGP_RUN_EXTRA_OPTIONS "-T sge -B '-V -b n -S /bin/bash' --pgp-batch-options-extra-prepare '-pe 1way 16 -q serial -l h_rt=01:30:00' --pgp-batch-options-extra-process '-pe 16way 64 -q normal -l h_rt=12:00:00'" CACHE STRING "Pass these optional arguments to pgp_run when testing. Tune this to your cluster configuration." FORCE)

set(MAKEFLOW_BUILD_MPI_BACKEND TRUE CACHE BOOL "Build MPI backend for Makeflow")
set(MAKEFLOW_MPI_PATH "$ENV{MPICH_HOME}" CACHE PATH "MPI Home dir for building Makeflow")
set(MAKEFLOW_GLOBUS_PATH "$ENV{GLOBUS_PATH}" CACHE PATH "GLOBUS Home dir for building Makeflow")

