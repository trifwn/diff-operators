function(define_targets)
    # -------------------------------------------------------------------------------------------------
    #                                            Operators
    # -------------------------------------------------------------------------------------------------
    add_library(operators_serial OBJECT ${SRC_OPERATORS}/operators_serial.f90)
    target_link_libraries(operators_serial PRIVATE )
    add_library(data_com OBJECT ${SRC_OPERATORS}/data_communication.f90)

    # -------------------------------------------------------------------------------------------------
    #                                            API Serial Vector Field Operators
    # -------------------------------------------------------------------------------------------------
    python_add_library(operators ${SRC_OPERATORS}/api_operators_serial.f90)
    target_link_libraries(operators PRIVATE operators_serial)
    set_target_properties(operators PROPERTIES
        OUTPUT_NAME "operators"
        # LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vpm_py_api
    )
    # -------------------------------------------------------------------------------------------------
    #                                           Operator Test Executable
    # -------------------------------------------------------------------------------------------------

    add_executable(test_operators ${SRC_TEST}/test_operators.f90)
    target_link_libraries(test_operators PUBLIC operators_serial data_com)

    # -------------------------------------------------------------------------------------------------
    #                                          Compiler Flags
    # -------------------------------------------------------------------------------------------------

    # Set Fortran compiler flags for all targets
    set_compiler_flags(data_com)
    set_compiler_flags(operators_serial)
    set_compiler_flags(test_operators)
endfunction()