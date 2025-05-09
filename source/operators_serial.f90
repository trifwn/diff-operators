module serial_vector_field_operators
    ! use omp_lib
    use iso_fortran_env, only: real64
    use custom_stencil_derivatives, only: compute_derivative
    implicit none

    
    private
    public :: divergence, curl, laplacian, gradient, calc_vector_laplacian_expansion
    public :: hessian, jacobian
    interface hessian
        module function calculate_scalar_hessian(f, dx, dy, dz) result(hes)
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in) :: dx, dy, dz
            real(real64), allocatable :: fx(:, :, :), fy(:, :, :), fz(:, :, :)
            real(real64), allocatable :: hes(:, :, :, :, :)
        end function calculate_scalar_hessian

        module function calculate_vector_hessian(f, dx, dy, dz) result(hes)
            real(real64), intent(in), target :: f(:, :, :, :)
            real(real64), intent(in) :: dx, dy, dz
            real(real64), allocatable :: fx(:, :, :), fy(:, :, :), fz(:, :, :)
            real(real64), allocatable :: hes(:, :, :, :, :, :)
        end function calculate_vector_hessian
    end interface hessian

    interface jacobian
        module function calculate_scalar_jacobian(f, dx, dy, dz) result(jac)
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in) :: dx, dy, dz
            real(real64), allocatable :: jac(:, :, :, :)
            integer :: nx, ny, nz
        end function calculate_scalar_jacobian

        module function calculate_vector_jacobian(f, dx, dy, dz) result(jac)
            real(real64), intent(in), target :: f(:, :, :, :)
            real(real64), intent(in) :: dx, dy, dz
            real(real64), allocatable :: jac(:, :, :, :, :)
            integer :: nx, ny, nz, i, neq
        end function calculate_vector_jacobian
    end interface jacobian

    interface laplacian
        module function calculate_scalar_laplacian(f, dx, dy, dz) result(result_lapl)
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in) :: dx, dy, dz
            real(real64), allocatable :: result_lapl(:, :, :)
            integer :: nx, ny, nz
        end function calculate_scalar_laplacian

        module function calculate_vector_laplacian(f, dx, dy, dz) result(result_lapl)
            real(real64), intent(in), target :: f(:, :, :, :)
            real(real64), intent(in) :: dx, dy, dz
            real(real64), allocatable :: result_lapl(:, :, :, :)
            integer :: neq, nx, ny, nz, i
        end function calculate_vector_laplacian
    end interface laplacian
contains

    module function calculate_scalar_jacobian(f, dx, dy, dz) result(jac)
        real(real64), intent(in), target :: f(:, :, :)
        real(real64), intent(in)         :: dx, dy, dz
        real(real64), allocatable :: jac(:, :, :, :)
        integer :: nx, ny, nz

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)

        allocate (jac(3, nx, ny, nz))
        jac(1, :, :, :) = compute_derivative(f, dx, dy, dz, 1, 1)
        jac(2, :, :, :) = compute_derivative(f, dx, dy, dz, 2, 1)
        jac(3, :, :, :) = compute_derivative(f, dx, dy, dz, 3, 1)
    end function calculate_scalar_jacobian

    module function calculate_vector_jacobian(f, dx, dy, dz) result(jac)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: jac(:, :, :, :, :)
        integer :: nx, ny, nz, i, neq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate (jac(neq, 3, nx, ny, nz))

        !$omp parallel do
        do i = 1, neq
            jac(i, :, :, :, :) = calculate_scalar_jacobian(f(i, :, :, :), dx, dy, dz)
        end do
        !$omp end parallel do
    end function calculate_vector_jacobian

    module function calculate_scalar_hessian(f, dx, dy, dz) result(hes)
        real(real64), intent(in), target :: f(:, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: hes(:, :, :, :, :)
        real(real64), allocatable :: fx(:, :, :), fy(:, :, :), fz(:, :, :)

        allocate (hes(3, 3, size(f, 1), size(f, 2), size(f, 3)))

        fx = compute_derivative(f, dx, dy, dz, 1, 1)
        fy = compute_derivative(f, dx, dy, dz, 2, 1)
        fz = compute_derivative(f, dx, dy, dz, 3, 1)

        ! Diagonal elements
        hes(1, 1, :, :, :) = compute_derivative(f, dx, dy, dz, 1, 2)
        hes(2, 2, :, :, :) = compute_derivative(f, dx, dy, dz, 2, 2)
        hes(3, 3, :, :, :) = compute_derivative(f, dx, dy, dz, 3, 2)

        ! Off-diagonal elements
        hes(1, 2, :, :, :) = compute_derivative(fx, dx, dy, dz, 2, 1)
        hes(1, 3, :, :, :) = compute_derivative(fx, dx, dy, dz, 3, 1)
        hes(2, 3, :, :, :) = compute_derivative(fy, dx, dy, dz, 3, 1)
        if (.true.) then ! Symmetric matrix
            hes(2, 1, :, :, :) = hes(1, 2, :, :, :)
            hes(3, 1, :, :, :) = hes(1, 3, :, :, :)
            hes(3, 2, :, :, :) = hes(2, 3, :, :, :)
        else             ! Non-symmetric matrix
            hes(2, 1, :, :, :) = compute_derivative(fy, dx, dy, dz, 1, 1)
            hes(3, 1, :, :, :) = compute_derivative(fz, dx, dy, dz, 1, 1)
            hes(3, 2, :, :, :) = compute_derivative(fz, dx, dy, dz, 2, 1)
        end if
        deallocate (fx, fy, fz)
    end function calculate_scalar_hessian

    module function calculate_vector_hessian(f, dx, dy, dz) result(hes)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: hes(:, :, :, :, :, :)
        integer :: nx, ny, nz, i, neq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate (hes(neq, 3, 3, nx, ny, nz))

        !$omp parallel do
        do i = 1, neq
            hes(i, :, :, :, :, :) = calculate_scalar_hessian(f(i, :, :, :), dx, dy, dz)
        end do
        !$omp end parallel do
    end function calculate_vector_hessian

    function divergence(f, dx, dy, dz) result(div_result)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: div_result(:, :, :)
        integer :: nx, ny, nz

        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        if (size(f, 1) /= 3) then
            print *, "Error: Divergence requires a 3D vector field"
            return
        end if

        allocate (div_result(nx, ny, nz))
        div_result(:, :, :) = compute_derivative(f(1, :, :, :), dx, dy, dz, 1, 1) + &
                              compute_derivative(f(2, :, :, :), dx, dy, dz, 2, 1) + &
                              compute_derivative(f(3, :, :, :), dx, dy, dz, 3, 1)
    end function divergence

    function curl(f, dx, dy, dz) result(result_curl)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: result_curl(:, :, :, :)
        real(real64), allocatable, target      :: J1(:, :, :), J2(:, :, :)
        integer :: nx, ny, nz

        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        if (size(f, 1) /= 3) then
            print *, "Error: Curl requires a 3D vector field"
            return
        end if

        allocate (result_curl(3, nx, ny, nz))
        ! We need to calculate the curl of a vector field F = (F1, F2, F3)
        ! That means calculating the offdiagnoal components of the Jacobian matrix
        ! of the vector field F. The Jacobian matrix is given by:
        ! J = | dF1/dx  dF1/dy  dF1/dz |
        !     | dF2/dx  dF2/dy  dF2/dz |
        !     | dF3/dx  dF3/dy  dF3/dz |

        ! The curl of a vector field is given by:
        ! curl(F) = (dF3/dy - dF2/dz, dF1/dz - dF3/dx, dF2/dx - dF1/dy)
        !         = (J32 - J23      , J13 - J31      , J21 - J12      )
        J1 = compute_derivative(f(3, :, :, :), dx, dy, dz, 2, 1) ! dF3/dy = J32
        J2 = compute_derivative(f(2, :, :, :), dx, dy, dz, 3, 1) ! dF2/dz = J23
        result_curl(1, :, :, :) = J1 - J2                          ! dF3/dy - dF2/dz = J32 - J23
        J1 = compute_derivative(f(1, :, :, :), dx, dy, dz, 3, 1) ! dF1/dz = J13
        J2 = compute_derivative(f(3, :, :, :), dx, dy, dz, 1, 1) ! dF3/dx = J31
        result_curl(2, :, :, :) = J1 - J2                          ! dF1/dz - dF3/dx = J13 - J31
        J1 = compute_derivative(f(2, :, :, :), dx, dy, dz, 1, 1) ! dF2/dx = J21
        J2 = compute_derivative(f(1, :, :, :), dx, dy, dz, 2, 1) ! dF1/dy = J12
        result_curl(3, :, :, :) = J1 - J2                          ! dF2/dx - dF1/dy = J21 - J12

        deallocate (J1, J2)
    end function curl

    function gradient(f, dx, dy, dz) result(result_grad)
        real(real64), intent(in), target :: f(:, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: result_grad(:, :, :, :)
        integer :: nx, ny, nz

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)

        allocate (result_grad(3, nx, ny, nz))
        result_grad(1, :, :, :) = compute_derivative(f(:, :, :), dx, dy, dz, 1, 1)
        result_grad(2, :, :, :) = compute_derivative(f(:, :, :), dx, dy, dz, 2, 1)
        result_grad(3, :, :, :) = compute_derivative(f(:, :, :), dx, dy, dz, 3, 1)
    end function gradient

    module function calculate_scalar_laplacian(f, dx, dy, dz) result(result_lapl)
        real(real64), intent(in), target :: f(:, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: result_lapl(:, :, :)
        integer :: nx, ny, nz

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)

        allocate (result_lapl(nx, ny, nz))
        result_lapl(:, :, :) = compute_derivative(f(:, :, :), dx, dy, dz, 1, 2) + &
                               compute_derivative(f(:, :, :), dx, dy, dz, 2, 2) + &
                               compute_derivative(f(:, :, :), dx, dy, dz, 3, 2)
    end function calculate_scalar_laplacian

    module function calculate_vector_laplacian(f, dx, dy, dz) result(result_lapl)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: result_lapl(:, :, :, :)
        integer :: neq, nx, ny, nz, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate (result_lapl(neq, nx, ny, nz))
        !$omp parallel do
        do i = 1, neq
            result_lapl(i, :, :, :) = compute_derivative(f(i, :, :, :), dx, dy, dz, 1, 2) + &
                                      compute_derivative(f(i, :, :, :), dx, dy, dz, 2, 2) + &
                                      compute_derivative(f(i, :, :, :), dx, dy, dz, 3, 2)
        end do
        !$omp end parallel do
    end function calculate_vector_laplacian

    function calc_vector_laplacian_expansion(f, dx, dy, dz) result(result_laplc)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), allocatable :: result_laplc(:, :, :, :), grad_div(:, :, :, :)
        real(real64), allocatable :: temp2(:, :, :), temp1(:, :, :, :)

        ! The correct Laplacian of a vector field F = (F1, F2, F3) is given by:
        ! Laplacian(F) = del(del · F) - del × (del × F)

        ! Calculate del · F
        temp2 = divergence(f, dx, dy, dz)

        ! Calculate del(del · F)
        grad_div = gradient(temp2, dx, dy, dz)

        ! Calculate del × F
        temp1 = curl(f, dx, dy, dz)

        ! Calculate del × (del × F)
        result_laplc = curl(temp1, dx, dy, dz)

        ! Combine the terms: del(del · F) - del × (del × F)
        result_laplc = grad_div - result_laplc
    end function calc_vector_laplacian_expansion
end module serial_vector_field_operators
