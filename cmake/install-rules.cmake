if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_INSTALL_INCLUDEDIR include/atomic-potentials CACHE PATH "")
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# find_package(<package>) call for consumers to find this project
set(package atomic-potentials)

install(
    TARGETS atomic-potentials_atomic-potentials
    EXPORT atomic-potentialsTargets
    RUNTIME COMPONENT atomic-potentials_Runtime
)

write_basic_package_version_file(
    "${package}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

# Allow package maintainers to freely override the path for the configs
set(
    atomic-potentials_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATADIR}/${package}"
    CACHE PATH "CMake package config location relative to the install prefix"
)
mark_as_advanced(atomic-potentials_INSTALL_CMAKEDIR)

install(
    FILES cmake/install-config.cmake
    DESTINATION "${atomic-potentials_INSTALL_CMAKEDIR}"
    RENAME "${package}Config.cmake"
    COMPONENT atomic-potentials_Development
)

install(
    FILES "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
    DESTINATION "${atomic-potentials_INSTALL_CMAKEDIR}"
    COMPONENT atomic-potentials_Development
)

install(
    EXPORT atomic-potentialsTargets
    NAMESPACE atomic-potentials::
    DESTINATION "${atomic-potentials_INSTALL_CMAKEDIR}"
    COMPONENT atomic-potentials_Development
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
