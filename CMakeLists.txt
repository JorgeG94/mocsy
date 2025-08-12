cmake_minimum_required(VERSION 3.22)
project(mocsy
  DESCRIPTION "Routines to model ocean carbonate system thermodynamics."
  LANGUAGES Fortran
  VERSION 2025.08.0)

include(CheckFortranCompilerFlag)
include(CMakePackageConfigHelpers)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Relwithdebinfo' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Relwithdebinfo"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
else()
  message(STATUS "Setting build type to ${CMAKE_BUILD_TYPE}")
endif()

set(MOCSY_PRECISION 2 CACHE STRING "Precision to use (1 or 2)")
option(MOCSY_PRECISION "Precision to use (1 or 2)" 2)
option(BUILD_SHARED_LIBS "Build shared/dynamic libraries" OFF)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

################################################################################

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -ffree-line-length-none -fno-range-check -Waliasing -Wampersand -Warray-bounds -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow -Wunused-parameter -Wintrinsic-shadow -Wno-align-commons -Wsurprising")
  set(FPE_EXCEPTIONS "-ffpe-trap=invalid,zero,overflow")

  check_fortran_compiler_flag("-fallow-invalid-boz" _boz_flag)
  check_fortran_compiler_flag("-fallow-argument-mismatch" _argmis_flag)
  if(_boz_flag)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-invalid-boz")
  endif()
  if(_argmis_flag)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -W -fbounds-check")

elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "INTEL" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -safe-cray-ptr -fpe0 -ftz -assume byterecl -i4 -traceback -nowarn -check noarg_temp_created -assume nobuffered_io -convert big_endian -grecord-gcc-switches -fp-model precise -fp-model source -align all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g3 -O0 -debug all -check -check noarg_temp_created -check nopointer -warn -warn noerrors -ftrapuv")

else()
  message(WARNING "Unsupported Fortran compiler with ID '${CMAKE_Fortran_COMPILER_ID}' will be used with default flags")
endif()

################################################################################
if(NOT TARGET test-drive::test-drive AND MOCSY_ENABLE_TESTING)
  enable_testing()
  find_package("test-drive" REQUIRED)
endif()

add_library(mocsy)

add_subdirectory(src)

set_target_properties(mocsy PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
  OUTPUT_NAME mocsy
  POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)

target_compile_definitions(mocsy PRIVATE USE_PRECISION=${MOCSY_PRECISION}) 

target_include_directories(mocsy PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS mocsy
  EXPORT mocsy-targets
  ARCHIVE DESTINATION lib
)

add_executable(test ${PROJECT_SOURCE_DIR}/test/test_mocsy.F90)
target_link_libraries(test PRIVATE mocsy)
target_include_directories(test PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:include>
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/mod/
  DESTINATION include
  FILES_MATCHING PATTERN "*.mod"
)

install(EXPORT mocsy-targets
  FILE mocsy-targets.cmake
  NAMESPACE mocsy::
  DESTINATION lib/cmake/mocsy
)

# Configure/install CMake package config and version file
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mocsy-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mocsy-config.cmake"
  INSTALL_DESTINATION lib/cmake/mocsy
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/mocsy-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/mocsy-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/mocsy-config-version.cmake"
  DESTINATION lib/cmake/mocsy
)

# Configure/install pkgconfig file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mocsy.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mocsy.pc"
  @ONLY
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/mocsy.pc"
  DESTINATION lib/pkgconfig
)

