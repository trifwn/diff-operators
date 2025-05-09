cmake_minimum_required(VERSION 3.22)
message(NOTICE "Configuring Operators API")

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# Define options
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add option for static build
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # -fPIC
option(BUILD_STATIC "Build static executable" ON)
option(CMAKE_INSTALL_LOCAL_ONLY "Install only the local project" ON)
set(CMAKE_INSTALL_LOCAL_ONLY ON)
SET(CMAKE_INSTALL_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/shared_libs")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Define the installation directories
if (CMAKE_INSTALL_PREFIX STREQUAL "" OR 
    NOT DEFINED CMAKE_INSTALL_PREFIX OR 
    CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT
)
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/build" CACHE PATH "install prefix" FORCE)
    message(NOTICE "SETTING CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endif()
message(NOTICE "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/../shared_libs")
set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/")

# Include custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FindCompilers)
find_compilers()

project(vpm VERSION 1.0)
      
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT "Fortran" IN_LIST languages)
    enable_language(Fortran)
endif()
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "Intel|GNU")
  message(WARNING "Unexpected Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# Setup compiler flags
include(CompilerFlags)
setup_compiler_flags()

# Setting directories
set(SRC_OPERATORS ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SRC_TEST ${CMAKE_CURRENT_SOURCE_DIR}/test)

# Determine the compiler name (e.g., "Intel", "GNU")
string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" COMPILER_NAME)
string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

set(LIB_DIRECTORY "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/${COMPILER_NAME}/lib/")
set(MOD_DIRECTORY "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/${COMPILER_NAME}/mod/")
set(BIN_DIRECTORY "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/${COMPILER_NAME}/bin/")

make_directory(${LIB_DIRECTORY})
make_directory(${MOD_DIRECTORY})
make_directory(${BIN_DIRECTORY})

# Set output directories
set(CMAKE_Fortran_MODULE_DIRECTORY ${MOD_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIRECTORY})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIRECTORY})


set(FISHPACK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../fishpack)

include(DefineTargets)
define_targets()

# -------------------------------------------------------------------------------------------------
#                                            Dependencies
# -------------------------------------------------------------------------------------------------

# Print the configuration summary 
message("-------------------------------------------------------------------------------------------------------------------")
message(NOTICE "Configuration summary:")
message(NOTICE "\tBUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(NOTICE "\tFortran_COMPILER: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")
message(NOTICE "\tFortran_COMPILER_VERSION: ${CMAKE_Fortran_COMPILER_VERSION}")
message(NOTICE "\tFortran_COMPILER_TARGET: ${CMAKE_Fortran_COMPILER_TARGET}")
message(NOTICE "\tFortran_COMPILER_EXTERNAL_TOOLCHAIN: ${CMAKE_Fortran_COMPILER_EXTERNAL_TOOLCHAIN}")
# Print the Fortran flags
message("Compiler flags:")
message(NOTICE "\tCMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(NOTICE "\tCMAKE_Fortran_FLAGS_DEBUG: ${Fortran_FLAGS_DEBUG}")
message(NOTICE "\tCMAKE_Fortran_FLAGS_RELEASE: ${Fortran_FLAGS_RELEASE}")
# Print the linker flags
message("Linker flags:")
message(NOTICE "\tCMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
message(NOTICE "\tCMAKE_EXE_LINKER_FLAGS_DEBUG: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message(NOTICE "\tCMAKE_EXE_LINKER_FLAGS_RELEASE: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
# Print the libraries
message("Output directories:")
message(NOTICE "\tCMAKE_Fortran_MODULE_DIRECTORY: ${CMAKE_Fortran_MODULE_DIRECTORY}")
message(NOTICE "\tCMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message(NOTICE "\tCMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(NOTICE "\tCMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Install the operators_api library
install(TARGETS operators
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)
