cmake_minimum_required(VERSION 3.13)

project("MyProject")

# Import HELYX targets
# NOTE: This assumes the environment has been sourced with activeBuild.shrc
# or similar - if not, an absolute path can be used
set(CMAKE_PREFIX_PATH "$ENV{HELYX_PROJECT_DIR}/platforms/$ENV{HELYX_OPTIONS}")
find_package(HELYX REQUIRED
OPTIONAL_COMPONENTS
HELYX-adjoint
HELYX-coupled
HELYX-hydro
HELYX-marine
# QUIET
CONFIG
)

# The following line supplies the HELYX CMake configuration
include(${HELYX_PROJECT_DIR}/etc/cmake/projectHeader.cmake)

set(TARGET_NAME "MyTarget-exe")

add_helyx_executable(${TARGET_NAME}
./app.C
)

helyx_additional_includes(${TARGET_NAME}
PRIVATE
surfMesh-includes
)

target_link_libraries(${TARGET_NAME}
PRIVATE finiteVolume
PRIVATE meshTools
PRIVATE OpenFOAM
)

# Set alternative output directory if desired
set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# The following line processes the targets in your project
include(${HELYX_PROJECT_DIR}/etc/cmake/projectFooter.cmake)
