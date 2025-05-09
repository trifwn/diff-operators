module custom_stencil_derivatives
    ! use omp_lib
    use iso_fortran_env, only: real64
    implicit none
    public :: compute_derivative

    interface compute_derivative
        module function derivative_1d(f, dx, dim, order, stencil) result(df)
            implicit none
            real(real64), intent(in), target :: f(:)
            real(real64), intent(in) :: dx
            integer, intent(in)          :: dim, order
            real(real64), optional, intent(in) :: stencil(:)
            real(real64), allocatable        :: df(:)
        end function

        module function vector_derivative_1d(f, dx, dim, order, stencil) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :)
            real(real64), intent(in)         :: dx
            integer, intent(in)          :: dim, order
            real(real64), optional, intent(in) :: stencil(:)
            real(real64), allocatable        :: df(:, :)
        end function

        module function derivative_2d(f, dx, dy, dim, order, stencil) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :)
            real(real64), intent(in)         :: dx, dy
            integer, intent(in)          :: dim, order
            real(real64), optional, intent(in) :: stencil(:)
            real(real64), allocatable        :: df(:, :)
        end function

        module function vector_derivative_2d(f, dx, dy, dim, order, stencil) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in)         :: dx, dy
            integer, intent(in)          :: dim, order
            real(real64), optional, intent(in) :: stencil(:)
            real(real64), allocatable        :: df(:, :, :)
        end function

        module function derivative_3d(f, dx, dy, dz, dim, order, stencil) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :, :)
            real(real64), intent(in)         :: dx, dy, dz
            integer, intent(in)          :: dim, order
            real(real64), optional, intent(in) :: stencil(:)
            real(real64), allocatable        :: df(:, :, :)
        end function

        module function vector_derivative_3d(f, dx, dy, dz, dim, order, stencil) result(df)
            implicit none
            real(real64), intent(in), target :: f(:, :, :, :)
            real(real64), intent(in)         :: dx, dy, dz
            integer, intent(in)          :: dim, order
            real(real64), optional, intent(in) :: stencil(:)
            real(real64), allocatable        :: df(:, :, :, :)
        end function
    end interface compute_derivative

    ! Private helper procedures
    private :: apply_stencil_1d
contains
    ! Helper function to apply arbitrary stencil to 1D data
    function apply_stencil_1d(f, dx, order, stencil, i) result(derivative)
        real(real64), intent(in) :: f(:)
        real(real64), intent(in) :: dx
        integer, intent(in)      :: order
        real(real64), intent(in) :: stencil(:)
        integer, intent(in)      :: i
        real(real64)             :: derivative
        
        integer :: half_width, j, idx
        real(real64) :: sum
        
        half_width = (size(stencil) - 1) / 2
        sum = 0.0_real64
        
        do j = -half_width, half_width
            idx = i + j
            ! Boundary handling (simple clamp to edge)
            if (idx < 1) idx = 1
            if (idx > size(f)) idx = size(f)
            
            sum = sum + stencil(j+half_width+1) * f(idx)
        end do
        
        derivative = sum / (dx**order)
    end function apply_stencil_1d

    module function derivative_1d(f, dx, dim, order, stencil) result(df)
        real(real64), intent(in), target :: f(:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: dim, order
        real(real64), optional, intent(in) :: stencil(:)
        real(real64), allocatable :: df(:)
        integer :: n, i
        real(real64), allocatable :: working_stencil(:)

        n = size(f)
        allocate(df(n))

        ! Use custom stencil if provided, otherwise use default stencils
        if (present(stencil)) then
            ! Apply custom stencil to each point
            do i = 1, n
                df(i) = apply_stencil_1d(f, dx, order, stencil, i)
            end do
        else
            ! Default stencils based on order
            select case (order)
            case (1)  ! First-order derivative (second-order accurate)
                ! Forward difference for the first point
                df(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*dx)
                
                ! Central difference for interior points (2nd order accurate)
                df(2:n-1) = (f(3:n) - f(1:n-2)) / (2*dx)
                
                ! Backward difference for the last point
                df(n) = (3*f(n) - 4*f(n-1) + f(n-2)) / (2*dx)
                
            case (2)  ! Second-order derivative
                ! Second-order accurate central difference for the first point
                df(1) = (f(3) - 2*f(2) + f(1))/(dx**2)
                
                ! Second-order accurate central difference for the interior points
                df(2:n-1) = (f(3:n) - 2*f(2:n-1) + f(1:n-2))/(dx**2)
                
                ! Second-order accurate central difference for the last point
                df(n) = (f(n) - 2*f(n-1) + f(n-2))/(dx**2)
            end select
        end if
    end function derivative_1d

    module function derivative_2d(f, dx, dy, dim, order, stencil) result(df)
        real(real64), intent(in), target :: f(:, :)
        real(real64), intent(in) :: dx, dy
        integer, intent(in) :: dim, order
        real(real64), optional, intent(in) :: stencil(:)
        real(real64), allocatable :: df(:, :)
        integer :: nx, ny, i, j

        nx = size(f, 1)
        ny = size(f, 2)
        allocate(df(nx, ny))

        !$omp parallel do 
        select case (dim)
        case (1)  ! x-derivative
            do j = 1, ny
                if (present(stencil)) then
                    df(:, j) = derivative_1d(f(:, j), dx, 1, order, stencil)
                else
                    df(:, j) = derivative_1d(f(:, j), dx, 1, order)
                end if
            end do
        case (2)  ! y-derivative
            do i = 1, nx
                if (present(stencil)) then
                    df(i, :) = derivative_1d(f(i, :), dy, 1, order, stencil)
                else
                    df(i, :) = derivative_1d(f(i, :), dy, 1, order)
                end if
            end do
        end select
        !$omp end parallel do
    end function derivative_2d

    module function derivative_3d(f, dx, dy, dz, dim, order, stencil) result(df)
        real(real64), intent(in), target  :: f(:, :, :)
        real(real64), intent(in) :: dx, dy, dz
        integer, intent(in) :: dim, order
        real(real64), optional, intent(in) :: stencil(:)
        real(real64), allocatable :: df(:, :, :)
        integer :: nx, ny, nz, i, j, k

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)

        allocate(df(nx, ny, nz))

        !$omp parallel do collapse(2)
        select case (dim)
        case (1)  ! x-derivative
            do k = 1, nz
                do j = 1, ny
                    if (present(stencil)) then
                        df(:, j, k) = derivative_1d(f(:, j, k), dx, 1, order, stencil)
                    else
                        df(:, j, k) = derivative_1d(f(:, j, k), dx, 1, order)
                    end if
                end do
            end do
        case (2)  ! y-derivative
            do k = 1, nz
                do i = 1, nx
                    if (present(stencil)) then
                        df(i, :, k) = derivative_1d(f(i, :, k), dy, 1, order, stencil)
                    else
                        df(i, :, k) = derivative_1d(f(i, :, k), dy, 1, order)
                    end if
                end do
            end do
        case (3)  ! z-derivative
            do j = 1, ny
                do i = 1, nx
                    if (present(stencil)) then
                        df(i, j, :) = derivative_1d(f(i, j, :), dz, 1, order, stencil)
                    else
                        df(i, j, :) = derivative_1d(f(i, j, :), dz, 1, order)
                    end if
                end do
            end do
        end select
        !$omp end parallel do
    end function derivative_3d

    module function vector_derivative_1d(f, dx, dim, order, stencil) result(df)
        real(real64), intent(in), target :: f(:, :)
        real(real64), intent(in) :: dx
        integer, intent(in) :: dim, order
        real(real64), optional, intent(in) :: stencil(:)
        real(real64), allocatable :: df(:, :)
        integer :: neq, n, i

        neq = size(f, 1)
        n = size(f, 2)
        allocate(df(neq, n))

        !$omp parallel do
        do i = 1, neq
            if (present(stencil)) then
                df(i, :) = derivative_1d(f(i, :), dx, dim, order, stencil)
            else
                df(i, :) = derivative_1d(f(i, :), dx, dim, order)
            end if
        end do
        !$omp end parallel do
    end function vector_derivative_1d

    module function vector_derivative_2d(f, dx, dy, dim, order, stencil) result(df)
        real(real64), intent(in), target :: f(:, :, :)
        real(real64), intent(in) ::  dx, dy
        integer, intent(in) :: dim, order
        real(real64), optional, intent(in) :: stencil(:)
        real(real64), allocatable :: df(:, :, :)
        integer :: neq, nx, ny, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        allocate(df(neq, nx, ny))

        !$omp parallel do
        do i = 1, neq
            select case (dim)
            case (1)  ! x-derivative
                if (present(stencil)) then
                    df(i, :, :) = derivative_2d(f(i, :, :), dx, dy, 1, order, stencil)
                else
                    df(i, :, :) = derivative_2d(f(i, :, :), dx, dy, 1, order)
                end if
            case (2)  ! y-derivative
                if (present(stencil)) then
                    df(i, :, :) = derivative_2d(f(i, :, :), dx, dy, 2, order, stencil)
                else
                    df(i, :, :) = derivative_2d(f(i, :, :), dx, dy, 2, order)
                end if
            end select
        end do
        !$omp end parallel do
    end function vector_derivative_2d

    module function vector_derivative_3d(f, dx, dy, dz, dim, order, stencil) result(df)
        real(real64), intent(in), target :: f(:, :, :, :)
        real(real64), intent(in) ::  dx, dy, dz
        integer, intent(in) :: dim, order
        real(real64), optional, intent(in) :: stencil(:)
        real(real64), allocatable :: df(:, :, :, :)
        integer :: neq, nx, ny, nz, eq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)
        allocate(df(neq, nx, ny, nz))

        !$omp parallel do
        do eq = 1, neq
            select case (dim)
            case (1)  ! x-derivative
                if (present(stencil)) then
                    df(eq, :, :, :) = derivative_3d(f(eq, :, :, :), dx, dy, dz, 1, order, stencil)
                else
                    df(eq, :, :, :) = derivative_3d(f(eq, :, :, :), dx, dy, dz, 1, order)
                end if
            case (2)  ! y-derivative
                if (present(stencil)) then
                    df(eq, :, :, :) = derivative_3d(f(eq, :, :, :), dx, dy, dz, 2, order, stencil)
                else
                    df(eq, :, :, :) = derivative_3d(f(eq, :, :, :), dx, dy, dz, 2, order)
                end if
            case (3)  ! z-derivative
                if (present(stencil)) then
                    df(eq, :, :, :) = derivative_3d(f(eq, :, :, :), dx, dy, dz, 3, order, stencil)
                else
                    df(eq, :, :, :) = derivative_3d(f(eq, :, :, :), dx, dy, dz, 3, order)
                end if
            end select
        end do
        !$omp end parallel do
    end function vector_derivative_3d
end module custom_stencil_derivatives