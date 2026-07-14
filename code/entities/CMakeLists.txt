cmake_minimum_required(VERSION 3.11)

project(revac_entities_portable C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

set(ENTITIES_DIR "${CMAKE_CURRENT_LIST_DIR}")
set(COMMON_DIR "${ENTITIES_DIR}/common")
set(CORE_DIR "${COMMON_DIR}/core")

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS
  "${CORE_DIR}/src/*.c")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/CMakeFiles/.*")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/build/.*")

add_library(revac_entities_core STATIC
  "${COMMON_DIR}/revocable_ac.c"
  "${ENTITIES_DIR}/ta/revac_ta.c"
  "${ENTITIES_DIR}/uav/revac_uav.c"
  "${ENTITIES_DIR}/gcs/revac_gcs.c"
  ${CORE_SOURCES})

target_include_directories(revac_entities_core PUBLIC
  "${COMMON_DIR}"
  "${ENTITIES_DIR}/ta"
  "${ENTITIES_DIR}/uav"
  "${ENTITIES_DIR}/gcs"
  "${CORE_DIR}/include"
  "${CORE_DIR}/include/arith"
  "${CORE_DIR}/include/arith/arith_q"
  "${CORE_DIR}/include/arith/arith_qiss"
  "${CORE_DIR}/include/arith/arith_qshow"
  "${CORE_DIR}/include/arith/arith_real"
  "${CORE_DIR}/include/arith/arith_z"
  "${CORE_DIR}/flint_headers")

target_compile_options(revac_entities_core PRIVATE
  -Wall -Wextra -Wundef -Wshadow -Wpointer-arith
  -Wcast-align -Wstrict-prototypes -Wwrite-strings
  -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code
  -march=native
  -O3)

target_link_libraries(revac_entities_core PUBLIC
  "${CORE_DIR}/lib/libflint.so"
  m)

set_target_properties(revac_entities_core PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(uav_demo
  "${ENTITIES_DIR}/uav/uav_demo.c")
target_link_libraries(uav_demo PRIVATE revac_entities_core)
set_target_properties(uav_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revocable_ac_demo
  "${ENTITIES_DIR}/test/revocable_ac_demo.c")
target_link_libraries(revocable_ac_demo PRIVATE revac_entities_core)
set_target_properties(revocable_ac_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

enable_testing()
add_test(NAME revocable_ac_demo COMMAND revocable_ac_demo)
