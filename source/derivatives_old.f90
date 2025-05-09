module derivatives
    ! use omp_lib
    use iso_fortran_env, only: real64
    implicit none
    public :: calculate_derivative

    interface calculate_derivative
        module function calculate_derivative_1d(f, dx, dim, order) result(df)
            implicit none
            real(real64), intent(in), target :: f(:)
            real(real64), intent(in) :: dx
            integer, intent(in)          :: dim, order
            real(real64), allocatable        :: df(:)
        end function

        module function calculate_vector_derivative_1d(f, dx, dim, order) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :)
            real(real64), intent(in)         :: dx
            integer, intent(in)          :: dim, order
            real(real64), allocatable        :: df(:, :)
        end function

        module function calculate_derivative_2d(f, dx, dy, dim, order) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :)
            real(real64), intent(in)         :: dx, dy
            integer, intent(in)          :: dim, order
            real(real64), allocatable        :: df(:, :)
        end function

        module function calculate_vector_derivative_2d(f, dx, dy, dim, order) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in)         :: dx, dy
            integer, intent(in)          :: dim, order
            real(real64), allocatable        :: df(:, :, :)
        end function

        module function calculate_derivative_3d(f, dx, dy, dz, dim, order) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in)         :: dx, dy, dz
            integer, intent(in)          :: dim, order
            real(real64), allocatable        :: df(:, :, :)
        end function

        module function calculate_vector_derivative_3d(f, dx, dy, dz, dim, order) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :, :, :)
            real(real64), intent(in)         :: dx, dy, dz
            integer, intent(in)          :: dim, order
            real(real64), allocatable        :: df(:, :, :, :)
        end function
    end interface calculate_derivative
contains
    module function calculate_derivative_1d(f, dx, dim, order) result(df)
        real(real64), intent(in), target :: f(:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: dim, order
        real(real64), allocatable :: df(:)
        integer :: n

        n = size(f)
        allocate (df(n))

        select case (order)
        case (1)  ! First-order derivative (second-order accurate)
            ! Forward difference for the first point
            df(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*dx)
            df(2) = (f(3) - f(1)) / (2*dx)

            ! Central difference for interior points (2nd order accurate)
            df(2:n - 1) = (f(3:n) - f(1:n - 2)) / (2*dx)
            ! Central difference for interior points (4th order accurate)
            ! df(3:n-2) = (f(1:n-4) - 8*f(2:n-3) + 8*f(4:n-1) - f(5:n)) / (12*dx)

            ! Backward difference for the last point
            df(n-1) = (f(n) - f(n - 2)) / (2*dx)
            df(n) = (3*f(n) - 4*f(n - 1) + f(n - 2)) / (2*dx)
        case (2)  ! Second-order derivative
            ! Second-order accurate central difference for the first and second points
            df(1) = (f(3) - 2*f(2) + f(1))/(dx**2)
            df(2) = (f(4) - 2*f(3) + f(2))/(dx**2)

            ! Fourth-order accurate central difference for the interior points
            ! df(3:n-2) = (-f(1:n-4) + 16*f(2:n-3) - 30*f(3:n-2) + 16*f(4:n-1) - f(5:n)) / (12*dx**2)
            ! Second-order accurate central difference for the interior points
            df(2:n - 1) = (f(3:n) - 2*f(2:n - 1) + f(1:n - 2))/(dx**2)

            ! Second-order accurate central difference for the last two points
            df(n - 1) = (f(n) - 2*f(n - 1) + f(n - 2))/(dx**2)
            df(n) = (f(n) - 2*f(n - 1) + f(n - 2))/(dx**2)
        end select

    end function calculate_derivative_1d

    module function calculate_derivative_2d(f, dx, dy, dim, order) result(df)
        real(real64), intent(in), target :: f(:, :)
        real(real64), intent(in) :: dx, dy
        integer, intent(in) :: dim, order
        real(real64), allocatable :: df(:, :)
        integer :: nx, ny, i, j

        nx = size(f, 1)
        ny = size(f, 2)
        allocate (df(nx, ny))

        !$omp parallel do 
        select case (dim)
        case (1)  ! x-derivative
            do j = 1, ny
                df(:, j) = calculate_derivative_1d(f(:, j), dx, 1, order)
            end do
        case (2)  ! y-derivative
            do i = 1, nx
                df(i, :) = calculate_derivative_1d(f(i, :), dy, 1, order)
            end do
        end select
        !$omp end parallel do
    end function calculate_derivative_2d

    module function calculate_derivative_3d(f, dx, dy, dz, dim, order) result(df)
        real(real64), intent(in), target  :: f(:, :, :)
        real(real64), intent(in) :: dx, dy, dz
        integer, intent(in) :: dim, order
        real(real64), allocatable :: df(:, :, :)
        integer :: nx, ny, nz, i, j, k

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)

        allocate (df(nx, ny, nz))

        !$omp parallel do collapse(2)
        select case (dim)
        case (1)  ! x-derivative
            do k = 1, nz
                do j = 1, ny
                    df(:, j, k) = calculate_derivative_1d(f(:, j, k), dx, 1, order)
                end do
            end do
        case (2)  ! y-derivative
            do k = 1, nz
                do i = 1, nx
                    df(i, :, k) = calculate_derivative_1d(f(i, :, k), dy, 1, order)
                end do
            end do
        case (3)  ! z-derivative
            do j = 1, ny
                do i = 1, nx
                    df(i, j, :) = calculate_derivative_1d(f(i, j, :), dz, 1, order)
                end do
            end do
        end select
        !$omp end parallel do
    end function calculate_derivative_3d

    module function calculate_vector_derivative_1d(f, dx, dim, order) result(df)
        real(real64), intent(in), target :: f(:, :)
        real(real64), intent(in) :: dx
        integer, intent(in) :: dim, order
        real(real64), allocatable :: df(:, :)
        integer :: neq, n, i

        neq = size(f, 1)
        n = size(f, 2)
        allocate (df(neq, n))

        !$omp parallel do
        do i = 1, neq
            df(i, :) = calculate_derivative_1d(f(i, :), dx, dim, order)
        end do
        !$omp end parallel do
    end function calculate_vector_derivative_1d

    module function calculate_vector_derivative_2d(f, dx, dy, dim, order) result(df)
        real(real64), intent(in), target :: f(:, :, :)
        real(real64), intent(in) ::  dx, dy
        integer, intent(in) :: dim, order
        real(real64), allocatable :: df(:, :, :)
        integer :: neq, nx, ny, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        allocate (df(neq, nx, ny))

        !$omp parallel do
        do i = 1, neq
            select case (dim)
            case (1)  ! x-derivative
                df(i, :, :) = calculate_derivative_2d(f(i, :, :), dx, dy, 1, order)
            case (2)  ! y-derivative
                df(i, :, :) = calculate_derivative_2d(f(i, :, :), dx, dy, 2, order)
            end select
        end do
        !$omp end parallel do
    end function calculate_vector_derivative_2d

    module function calculate_vector_derivative_3d(f, dx, dy, dz, dim, order) result(df)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) ::  dx, dy, dz
        integer, intent(in) :: dim, order
        real(real64), allocatable :: df(:, :, :, :)
        integer :: neq, nx, ny, nz, eq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)
        allocate (df(neq, nx, ny, nz))

        !$omp parallel do
        do eq = 1, neq
            select case (dim)
            case (1)  ! x-derivative
                df(eq, :, :, :) = calculate_derivative_3d(f(eq, :, :, :), dx, dy, dz, 1, order)
            case (2)  ! y-derivative
                df(eq, :, :, :) = calculate_derivative_3d(f(eq, :, :, :), dx, dy, dz, 2, order)
            case (3)  ! z-derivative
                df(eq, :, :, :) = calculate_derivative_3d(f(eq, :, :, :), dx, dy, dz, 3, order)
            end select
        end do
        !$omp end parallel do

    end function calculate_vector_derivative_3d
end module derivatives