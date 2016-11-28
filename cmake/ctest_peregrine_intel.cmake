# -----------------------------------------------------------
# -- Variables for intel ctesting on Peregrine
# -----------------------------------------------------------

set(compiler                            "intel")
set(CONFIG_FILE                         "do-configNalu_release-${compiler}")

# -----------------------------------------------------------
# -- Include common test config
# -----------------------------------------------------------
include(../cmake/ctest_common.cmake)
