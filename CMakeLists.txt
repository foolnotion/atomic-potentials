cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    atomic-potentials
    VERSION 0.1.0
    DESCRIPTION "Symbolic regression modeling of interatomic potentials"
    HOMEPAGE_URL "https://gitlab.fh-ooe.at/jrz-symreg/atomic-potentials"
    LANGUAGES CXX
)

docs_early_return()

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    atomic-potentials_lib OBJECT
    source/lib.cpp
)

target_include_directories(
    atomic-potentials_lib ${atomic-potentials_warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(atomic-potentials_lib PUBLIC cxx_std_17)

# ---- Declare executable ----

add_executable(atomic-potentials_atomic-potentials source/main.cpp)
add_executable(atomic-potentials::atomic-potentials ALIAS atomic-potentials_atomic-potentials)

set_target_properties(
    atomic-potentials_atomic-potentials PROPERTIES
    OUTPUT_NAME atomic-potentials
    EXPORT_NAME atomic-potentials
)

target_compile_features(atomic-potentials_atomic-potentials PRIVATE cxx_std_17)

target_link_libraries(atomic-potentials_atomic-potentials PRIVATE atomic-potentials_lib)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT atomic-potentials_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of atomic-potentials"
  )
endif()

include(cmake/dev-mode.cmake)
