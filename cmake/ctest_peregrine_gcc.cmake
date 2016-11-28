# -----------------------------------------------------------
# -- Variables for gcc ctesting on Peregrine
# -----------------------------------------------------------

set(compiler                            "gcc")
set(CONFIG_FILE                         "do-configNalu_release-${compiler}")

# -----------------------------------------------------------
# -- Include common test config
# -----------------------------------------------------------
include(../cmake/ctest_common.cmake)
