cmake_minimum_required(VERSION 3.16)
project(NEP_CPU LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(NEP_CPU_OPENMP "Enable OpenMP" OFF)
option(NEP_CPU_USE_TABLE "Use tabulated radial functions" OFF)

add_library(nep_cpu
  src/nep.cpp
  src/nep_model_file.cpp
  src/ewald_nep.cpp
  src/neighbor_nep.cpp
)
add_library(nep::nep_cpu ALIAS nep_cpu)
target_compile_features(nep_cpu PUBLIC cxx_std_11)
target_include_directories(nep_cpu PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/NEP_CPU>
)
set_target_properties(nep_cpu PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)

if(NEP_CPU_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(nep_cpu PUBLIC OpenMP::OpenMP_CXX)
endif()

if(NEP_CPU_USE_TABLE)
  target_compile_definitions(nep_cpu PUBLIC USE_TABLE_FOR_RADIAL_FUNCTIONS)
endif()

install(TARGETS nep_cpu EXPORT NEP_CPUTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES src/nep.h src/ewald_nep.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/NEP_CPU
)
install(EXPORT NEP_CPUTargets
  NAMESPACE nep::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NEP_CPU
)
configure_package_config_file(cmake/NEP_CPUConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/NEP_CPUConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NEP_CPU
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NEP_CPUConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NEP_CPU
)
