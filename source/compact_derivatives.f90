module compact_derivatives
    use iso_fortran_env, only: real64
    implicit none

    private
    public :: compute_compact_derivative

    ! Interfaces for the main compact derivative computation functions
    interface compute_compact_derivative
        module function compute_compact_derivative_1d(f, dx, order, lhs_stencil, rhs_stencil, periodic) result(df)
            implicit none
            real(real64), intent(in) :: f(:)
            real(real64), intent(in) :: dx
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)  ! [alpha_{i-1}, alpha_i, alpha_{i+1}]
            real(real64), intent(in) :: rhs_stencil(:)  ! [a_i-2, a_i-1, a_i, a_i+1, a_i+2]
            logical, intent(in), optional :: periodic
            real(real64), allocatable :: df(:)
        end function

        module function compute_compact_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:)
            real(real64), intent(in) :: dx, dy
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            logical, intent(in), optional :: periodic
            real(real64), allocatable :: df(:,:)
        end function

        module function compute_compact_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:,:)
            real(real64), intent(in) :: dx, dy, dz
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            logical, intent(in), optional :: periodic
            real(real64), allocatable :: df(:,:,:)
        end function

        module function compute_compact_vector_derivative_1d(f, dx, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:)
            real(real64), intent(in) :: dx
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            logical, intent(in), optional :: periodic
            real(real64), allocatable :: df(:,:)
        end function

        module function compute_compact_vector_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:,:)
            real(real64), intent(in) :: dx, dy
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            logical, intent(in), optional :: periodic
            real(real64), allocatable :: df(:,:,:)
        end function

        module function compute_compact_vector_derivative_3d( &
                f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:,:,:)
            real(real64), intent(in) :: dx, dy, dz
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            logical, intent(in), optional :: periodic
            real(real64), allocatable :: df(:,:,:,:)
        end function
    end interface compute_compact_derivative

    ! Private helper functions
    private :: tdma_solver, cyclic_tdma_solver

contains
    ! Thomas algorithm (TDMA - Tridiagonal Matrix Algorithm) solver
    ! Solves system Ax=d where A is tridiagonal with diagonals a, b, c
    subroutine tdma_solver(a, b, c, d, x, n)
        integer, intent(in) :: n
        real(real64), intent(in) :: a(n), b(n), c(n), d(n)
        real(real64), intent(out) :: x(n)

        real(real64) :: cp(n), dp(n)
        integer :: i

        ! Forward sweep - elimination
        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)

        do i = 2, n
            cp(i) = c(i) / (b(i) - a(i) * cp(i-1))
            dp(i) = (d(i) - a(i) * dp(i-1)) / (b(i) - a(i) * cp(i-1))
        end do

        ! Backward substitution
        x(n) = dp(n)
        do i = n-1, 1, -1
            x(i) = dp(i) - cp(i) * x(i+1)
        end do
    end subroutine tdma_solver

    ! Cyclic TDMA solver using Sherman-Morrison formula
    ! Solves a cyclic tridiagonal system where alpha = A(1,n) and beta = A(n,1)
    subroutine cyclic_tdma_solver(a, b, c, d, x, n, alpha, beta)
        integer, intent(in) :: n
        real(real64), intent(in) :: a(n), b(n), c(n), d(n)
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(out) :: x(n)

        real(real64) :: gamma, factor
        real(real64) :: bb(n), u(n), y(n), z(n)
        integer :: i

        gamma = -b(1)

        ! Modified diagonal for the non-cyclic sub-problem
        bb = b
        bb(1) = b(1) - gamma
        bb(n) = b(n) - alpha * beta / gamma

        ! Solve A_tri * y = d
        call tdma_solver(a, bb, c, d, y, n)

        ! Solve A_tri * z = u where u = [gamma, 0, ..., 0, alpha]
        u = 0.0_real64
        u(1) = gamma
        u(n) = alpha
        call tdma_solver(a, bb, c, u, z, n)

        ! Apply Sherman-Morrison correction
        ! v = [1, 0, ..., 0, beta/gamma]
        factor = (y(1) + (beta / gamma) * y(n)) / (1.0_real64 + z(1) + (beta / gamma) * z(n))

        do i = 1, n
            x(i) = y(i) - factor * z(i)
        end do
    end subroutine cyclic_tdma_solver

    ! Helper function to apply the RHS stencil for a specific point
    function apply_rhs_stencil(f, i, n, rhs_stencil, dx, order, is_periodic) result(rhs_val)
        real(real64), intent(in) :: f(:)
        integer, intent(in) :: i, n
        real(real64), intent(in) :: rhs_stencil(:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: order
        logical, intent(in) :: is_periodic
        real(real64) :: rhs_val

        integer :: j, idx, stencil_size, offset

        stencil_size = size(rhs_stencil)
        offset = stencil_size / 2

        rhs_val = 0.0_real64

        do j = 1, stencil_size
            idx = i + (j - offset - 1)

            if (is_periodic) then
                ! Periodic wrapping using modulo (handles negative indices correctly)
                idx = modulo(idx - 1, n) + 1
            else
                ! Clamp to boundary for near-boundary interior points
                if (idx < 1) idx = 1
                if (idx > n) idx = n
            end if

            rhs_val = rhs_val + rhs_stencil(j) * f(idx)
        end do

        ! Scale by appropriate power of dx
        rhs_val = rhs_val / (dx**order)
    end function apply_rhs_stencil

    ! 1D compact derivative computation
    module function compute_compact_derivative_1d(f, dx, order, lhs_stencil, rhs_stencil, periodic) result(df)
        real(real64), intent(in) :: f(:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)  ! [alpha_{i-1}, alpha_i, alpha_{i+1}]
        real(real64), intent(in) :: rhs_stencil(:)  ! [a_i-2, a_i-1, a_i, a_i+1, a_i+2]
        logical, intent(in), optional :: periodic
        real(real64), allocatable :: df(:)

        integer :: n, i
        real(real64), allocatable :: a(:), b(:), c(:), d(:)
        logical :: is_periodic

        is_periodic = .false.
        if (present(periodic)) is_periodic = periodic

        n = size(f)
        allocate(df(n))
        allocate(a(n), b(n), c(n), d(n))

        ! Set up tridiagonal system coefficients for all points
        do i = 1, n
            a(i) = lhs_stencil(1)
            b(i) = lhs_stencil(2)
            c(i) = lhs_stencil(3)
            d(i) = apply_rhs_stencil(f, i, n, rhs_stencil, dx, order, is_periodic)
        end do

        if (is_periodic) then
            ! Periodic domain: solve cyclic tridiagonal system
            call cyclic_tdma_solver(a, b, c, d, df, n, lhs_stencil(1), lhs_stencil(3))
        else
            ! Non-periodic domain: apply proper boundary closures
            if (order == 1 .and. n >= 3) then
                ! Lele's 3rd-order boundary closure for first derivative
                ! Left boundary: f'_1 + 2*f'_2 = (-5*f_1 + 4*f_2 + f_3) / (2*dx)
                a(1) = 0.0_real64
                b(1) = 1.0_real64
                c(1) = 2.0_real64
                d(1) = (-5.0_real64*f(1) + 4.0_real64*f(2) + f(3)) / (2.0_real64*dx)

                ! Right boundary: 2*f'_{n-1} + f'_n = (-f_{n-2} - 4*f_{n-1} + 5*f_n) / (2*dx)
                a(n) = 2.0_real64
                b(n) = 1.0_real64
                c(n) = 0.0_real64
                d(n) = (-f(n-2) - 4.0_real64*f(n-1) + 5.0_real64*f(n)) / (2.0_real64*dx)
            else if (order == 2 .and. n >= 4) then
                ! Explicit 2nd-order one-sided formulas for second derivative
                ! Left boundary: f''_1 = (2*f_1 - 5*f_2 + 4*f_3 - f_4) / dx^2
                a(1) = 0.0_real64
                b(1) = 1.0_real64
                c(1) = 0.0_real64
                d(1) = (2.0_real64*f(1) - 5.0_real64*f(2) + 4.0_real64*f(3) - f(4)) / (dx**2)

                ! Right boundary: f''_n = (-f_{n-3} + 4*f_{n-2} - 5*f_{n-1} + 2*f_n) / dx^2
                a(n) = 0.0_real64
                b(n) = 1.0_real64
                c(n) = 0.0_real64
                d(n) = (-f(n-3) + 4.0_real64*f(n-2) - 5.0_real64*f(n-1) + 2.0_real64*f(n)) / (dx**2)
            else
                ! Fallback: decouple boundary points from the system
                a(1) = 0.0_real64
                c(n) = 0.0_real64
            end if

            call tdma_solver(a, b, c, d, df, n)
        end if
    end function compute_compact_derivative_1d

    ! 2D compact derivative computation
    module function compute_compact_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
        real(real64), intent(in) :: f(:,:)
        real(real64), intent(in) :: dx, dy
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        logical, intent(in), optional :: periodic
        real(real64), allocatable :: df(:,:)

        integer :: nx, ny, i, j
        real(real64), allocatable :: temp(:)

        nx = size(f, 1)
        ny = size(f, 2)
        allocate(df(nx, ny))

        select case(dim)
        case(1)  ! x-derivative
            !$omp parallel private(j, temp)
            allocate(temp(nx))
            !$omp do
            do j = 1, ny
                temp = compute_compact_derivative_1d(f(:,j), dx, order, lhs_stencil, rhs_stencil, periodic)
                df(:,j) = temp
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel

        case(2)  ! y-derivative
            !$omp parallel private(i, temp)
            allocate(temp(ny))
            !$omp do
            do i = 1, nx
                temp = compute_compact_derivative_1d(f(i,:), dy, order, lhs_stencil, rhs_stencil, periodic)
                df(i,:) = temp
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel
        end select
    end function compute_compact_derivative_2d

    ! 3D compact derivative computation
    module function compute_compact_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
        real(real64), intent(in) :: f(:,:,:)
        real(real64), intent(in) :: dx, dy, dz
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        logical, intent(in), optional :: periodic
        real(real64), allocatable :: df(:,:,:)

        integer :: nx, ny, nz, i, j, k
        real(real64), allocatable :: temp(:)

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)
        allocate(df(nx, ny, nz))

        select case(dim)
        case(1)  ! x-derivative
            !$omp parallel private(j, k, temp)
            allocate(temp(nx))
            !$omp do collapse(2)
            do k = 1, nz
                do j = 1, ny
                    temp = compute_compact_derivative_1d(f(:,j,k), dx, order, lhs_stencil, rhs_stencil, periodic)
                    df(:,j,k) = temp
                end do
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel

        case(2)  ! y-derivative
            !$omp parallel private(i, k, temp)
            allocate(temp(ny))
            !$omp do collapse(2)
            do k = 1, nz
                do i = 1, nx
                    temp = compute_compact_derivative_1d(f(i,:,k), dy, order, lhs_stencil, rhs_stencil, periodic)
                    df(i,:,k) = temp
                end do
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel

        case(3)  ! z-derivative
            !$omp parallel private(i, j, temp)
            allocate(temp(nz))
            !$omp do collapse(2)
            do j = 1, ny
                do i = 1, nx
                    temp = compute_compact_derivative_1d(f(i,j,:), dz, order, lhs_stencil, rhs_stencil, periodic)
                    df(i,j,:) = temp
                end do
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel
        end select
    end function compute_compact_derivative_3d

    ! Vector version of 1D compact derivative
    module function compute_compact_vector_derivative_1d(f, dx, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
        real(real64), intent(in) :: f(:,:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        logical, intent(in), optional :: periodic
        real(real64), allocatable :: df(:,:)

        integer :: neq, n, i

        neq = size(f, 1)
        n = size(f, 2)
        allocate(df(neq, n))

        !$omp parallel do
        do i = 1, neq
            df(i,:) = compute_compact_derivative_1d(f(i,:), dx, order, lhs_stencil, rhs_stencil, periodic)
        end do
        !$omp end parallel do
    end function compute_compact_vector_derivative_1d

    ! Vector version of 2D compact derivative
    module function compute_compact_vector_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
        real(real64), intent(in) :: f(:,:,:)
        real(real64), intent(in) :: dx, dy
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        logical, intent(in), optional :: periodic
        real(real64), allocatable :: df(:,:,:)

        integer :: neq, nx, ny, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        allocate(df(neq, nx, ny))

        !$omp parallel do
        do i = 1, neq
            df(i,:,:) = compute_compact_derivative_2d(f(i,:,:), dx, dy, dim, order, lhs_stencil, rhs_stencil, periodic)
        end do
        !$omp end parallel do
    end function compute_compact_vector_derivative_2d

    ! Vector version of 3D compact derivative
    module function compute_compact_vector_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil, periodic) result(df)
        real(real64), intent(in) :: f(:,:,:,:)
        real(real64), intent(in) :: dx, dy, dz
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        logical, intent(in), optional :: periodic
        real(real64), allocatable :: df(:,:,:,:)

        integer :: neq, nx, ny, nz, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)
        allocate(df(neq, nx, ny, nz))

        !$omp parallel do
        do i = 1, neq
            df(i,:,:,:) = compute_compact_derivative_3d(f(i,:,:,:), dx, dy, dz, dim, order, lhs_stencil, rhs_stencil, periodic)
        end do
        !$omp end parallel do
    end function compute_compact_vector_derivative_3d

end module compact_derivatives