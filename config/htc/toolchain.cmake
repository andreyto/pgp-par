##CMake is quite difficult (if not impossible) when you need to pass a shell-quoted string
##like below to some subprocess unchanged. The quotations and STRING variable
##type below do work (with extra shlex splitting done in the test Python code).
##Follow this example when editing the value, or you will break the test.
#set(PGP_TEST_PGP_RUN_EXTRA_OPTIONS "-T sge -B '-P 0534 -l fast -b n -S /bin/bash'" CACHE STRING "Pass these optional arguments to pgp_run when testing. Tune this to your cluster configuration." FORCE)


