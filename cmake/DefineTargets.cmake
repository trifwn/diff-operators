function(define_targets)
    # -------------------------------------------------------------------------------------------------
    #                                            Operators
    # -------------------------------------------------------------------------------------------------
    add_library(derivatives_old OBJECT ${SRC_OPERATORS}/derivatives_old.f90)

    add_library(custom_stencil_derivatives OBJECT ${SRC_OPERATORS}/custom_stencil_derivatives.f90)
    add_library(compact_derivatives OBJECT ${SRC_OPERATORS}/compact_derivatives.f90)

    add_library(operators_serial OBJECT ${SRC_OPERATORS}/operators_serial.f90)
    target_link_libraries(operators_serial PRIVATE custom_stencil_derivatives)
    add_library(data_com OBJECT ${SRC_OPERATORS}/data_communication.f90)

    # -------------------------------------------------------------------------------------------------
    #                                            API Serial Vector Field Operators
    # -------------------------------------------------------------------------------------------------
    python_add_library(operators ${SRC_OPERATORS}/api_operators_serial.f90)
    target_link_libraries(operators PRIVATE operators_serial custom_stencil_derivatives compact_derivatives)
    set_target_properties(operators PROPERTIES
        OUTPUT_NAME "operators"
    )
    # -------------------------------------------------------------------------------------------------
    #                                           Operator Test Executable
    # -------------------------------------------------------------------------------------------------

    add_executable(test_operators ${SRC_TEST}/test_operators.f90)
    target_link_libraries(test_operators PUBLIC custom_stencil_derivatives operators_serial data_com)


    add_executable(test_compact ${SRC_TEST}/test_compact.f90)
    target_link_libraries(test_compact PUBLIC compact_derivatives operators_serial custom_stencil_derivatives)

    # -------------------------------------------------------------------------------------------------
    #                                          Compiler Flags
    # -------------------------------------------------------------------------------------------------

    # Set Fortran compiler flags for all targets
    set_compiler_flags(data_com)
    set_compiler_flags(operators_serial)
    set_compiler_flags(test_operators)
endfunction()