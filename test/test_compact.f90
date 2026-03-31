program test_compact_derivatives
    use iso_fortran_env, only: real64, output_unit
    use compact_derivatives
    use custom_stencil_derivatives, only: compute_derivative  ! Optional: for comparison
    implicit none

    ! Parameters for testing
    integer, parameter :: n_points = 100
    integer :: i, test_case
    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    real(real64) :: domain_length, dx, dx_np
    integer, parameter :: n_np = 101  ! Non-periodic grid size

    ! Arrays for grid and function values
    real(real64), allocatable :: x(:), f(:), df_exact(:)
    real(real64), allocatable :: df_compact(:), df_fd(:)
    ! Non-periodic arrays
    real(real64), allocatable :: x_np(:), f_np(:), df_exact_np(:)
    real(real64), allocatable :: df_compact_np(:), df_fd_np(:)

    ! Various compact scheme stencils
    real(real64) :: lhs_4th_order(3), rhs_4th_order(5)
    real(real64) :: lhs_6th_order(3), rhs_6th_order(5)
    real(real64) :: lhs_2nd_deriv(3), rhs_2nd_deriv(5)

    ! Error metrics
    real(real64) :: l2_error_compact, l2_error_fd
    real(real64) :: max_error_compact, max_error_fd

    ! File name for output
    character(len=100) :: filename

    ! Initialize stencils for compact schemes

    ! 4th-order compact scheme for first derivative (Pade, alpha=1/4)
    lhs_4th_order = [0.25_real64, 1.0_real64, 0.25_real64]
    rhs_4th_order = [0.0_real64, -0.75_real64, 0.0_real64, 0.75_real64, 0.0_real64]

    ! 6th-order compact scheme for first derivative (alpha=1/3)
    lhs_6th_order = [1.0_real64/3, 1.0_real64, 1.0_real64/3]
    rhs_6th_order = [-1.0_real64/36, -7.0_real64/9, 0.0_real64, 7.0_real64/9, 1.0_real64/36]

    ! 4th-order compact scheme for second derivative (alpha=1/10, a=6/5)
    lhs_2nd_deriv = [0.1_real64, 1.0_real64, 0.1_real64]
    rhs_2nd_deriv = [0.0_real64, 1.2_real64, -2.4_real64, 1.2_real64, 0.0_real64]

    ! Setup periodic grid (excludes duplicate endpoint)
    domain_length = 2.0_real64 * pi
    dx = domain_length / n_points

    allocate(x(n_points), f(n_points), df_exact(n_points))
    allocate(df_compact(n_points), df_fd(n_points))

    do i = 1, n_points
        x(i) = (i - 1) * dx
    end do

    ! Setup non-periodic grid (includes both endpoints)
    dx_np = domain_length / (n_np - 1)

    allocate(x_np(n_np), f_np(n_np), df_exact_np(n_np))
    allocate(df_compact_np(n_np), df_fd_np(n_np))

    do i = 1, n_np
        x_np(i) = (i - 1) * dx_np
    end do

    ! Run different test cases
    do test_case = 1, 4
        ! Setup test function and its exact derivative based on test case
        select case(test_case)
        case(1)
            print *, "Test Case 1: f(x) = sin(x), first derivative (periodic)"
            f = sin(x)
            df_exact = cos(x)

            ! Compute first derivative using compact scheme (4th order, periodic)
            df_compact = compute_compact_derivative(f, dx, 1, lhs_4th_order, rhs_4th_order, periodic=.true.)

            ! Compare with standard finite difference
            df_fd = compute_derivative(f, dx, 1, 1)

            ! Compute error metrics
            l2_error_compact = sqrt(sum((df_compact - df_exact)**2) / n_points)
            max_error_compact = maxval(abs(df_compact - df_exact))

            l2_error_fd = sqrt(sum((df_fd - df_exact)**2) / n_points)
            max_error_fd = maxval(abs(df_fd - df_exact))

            ! Print results
            call print_results(n_points, dx, l2_error_compact, max_error_compact, l2_error_fd, max_error_fd)

            ! Output data to file
            write(filename, "(A,I3.3,A)") "test_case", test_case, ".dat"
            call write_results(filename, x, f, df_exact, df_compact, df_fd)

        case(2)
            print *, "Test Case 2: f(x) = sin(x), second derivative (periodic)"
            f = sin(x)
            df_exact = -sin(x)

            ! Compute second derivative using compact scheme (periodic)
            df_compact = compute_compact_derivative(f, dx, 2, lhs_2nd_deriv, rhs_2nd_deriv, periodic=.true.)

            ! Compare with standard finite difference
            df_fd = compute_derivative(f, dx, 1, 2)

            ! Compute error metrics
            l2_error_compact = sqrt(sum((df_compact - df_exact)**2) / n_points)
            max_error_compact = maxval(abs(df_compact - df_exact))

            l2_error_fd = sqrt(sum((df_fd - df_exact)**2) / n_points)
            max_error_fd = maxval(abs(df_fd - df_exact))

            call print_results(n_points, dx, l2_error_compact, max_error_compact, l2_error_fd, max_error_fd)

            write(filename, "(A,I3.3,A)") "test_case", test_case, ".dat"
            call write_results(filename, x, f, df_exact, df_compact, df_fd)

        case(3)
            print *, "Test Case 3: f(x) = exp(sin(x)), first derivative (periodic)"
            f = exp(sin(x))
            df_exact = cos(x) * exp(sin(x))

            ! Compute first derivative using 6th-order compact scheme (periodic)
            df_compact = compute_compact_derivative(f, dx, 1, lhs_6th_order, rhs_6th_order, periodic=.true.)

            ! Compare with standard finite difference
            df_fd = compute_derivative(f, dx, 1, 1)

            ! Compute error metrics
            l2_error_compact = sqrt(sum((df_compact - df_exact)**2) / n_points)
            max_error_compact = maxval(abs(df_compact - df_exact))

            l2_error_fd = sqrt(sum((df_fd - df_exact)**2) / n_points)
            max_error_fd = maxval(abs(df_fd - df_exact))

            call print_results(n_points, dx, l2_error_compact, max_error_compact, l2_error_fd, max_error_fd)

            write(filename, "(A,I3.3,A)") "test_case", test_case, ".dat"
            call write_results(filename, x, f, df_exact, df_compact, df_fd)

        case(4)
            print *, "Test Case 4: f(x) = x^3, first derivative (non-periodic)"
            f_np = x_np**3
            df_exact_np = 3.0_real64 * x_np**2

            ! Compute first derivative using 4th-order compact scheme (non-periodic)
            df_compact_np = compute_compact_derivative(f_np, dx_np, 1, lhs_4th_order, rhs_4th_order)

            ! Compare with standard finite difference
            df_fd_np = compute_derivative(f_np, dx_np, 1, 1)

            ! Compute error metrics (excluding boundary points for fairness)
            l2_error_compact = sqrt(sum((df_compact_np(3:n_np-2) - df_exact_np(3:n_np-2))**2) / (n_np-4))
            max_error_compact = maxval(abs(df_compact_np(3:n_np-2) - df_exact_np(3:n_np-2)))

            l2_error_fd = sqrt(sum((df_fd_np(3:n_np-2) - df_exact_np(3:n_np-2))**2) / (n_np-4))
            max_error_fd = maxval(abs(df_fd_np(3:n_np-2) - df_exact_np(3:n_np-2)))

            call print_results(n_np, dx_np, l2_error_compact, max_error_compact, l2_error_fd, max_error_fd)

            write(filename, "(A,I3.3,A)") "test_case", test_case, ".dat"
            call write_results(filename, x_np, f_np, df_exact_np, df_compact_np, df_fd_np)
        end select

        print *, ""
    end do

    ! Test 2D derivatives
    call test_2d_derivative()

    deallocate(x, f, df_exact, df_compact, df_fd)
    deallocate(x_np, f_np, df_exact_np, df_compact_np, df_fd_np)

contains
    ! Subroutine to print error results
    subroutine print_results(npts, h, l2_compact, max_compact, l2_fd, max_fd)
        integer, intent(in) :: npts
        real(real64), intent(in) :: h, l2_compact, max_compact, l2_fd, max_fd

        print *, "  Grid points: ", npts
        print *, "  Grid spacing: ", h
        print *, "  Compact scheme results:"
        print *, "    L2 error: ", l2_compact
        print *, "    Max error: ", max_compact
        print *, "  Finite difference results:"
        print *, "    L2 error: ", l2_fd
        print *, "    Max error: ", max_fd
        if (l2_compact > 0.0_real64) then
            print *, "  Error ratio (FD/Compact): ", l2_fd / l2_compact
        end if
    end subroutine print_results

    ! Subroutine to test 2D derivatives
    subroutine test_2d_derivative()
        integer, parameter :: nx = 50, ny = 50
        real(real64) :: domain_x, domain_y
        real(real64) :: dx2, dy2
        real(real64), allocatable :: x2d(:), y2d(:), f2d(:,:), df2d_exact(:,:), df2d_compact(:,:)
        integer :: i, j
        real(real64) :: l2_error, max_error

        print *, "Test Case 5: 2D derivative - f(x,y) = sin(x)*cos(y) (periodic)"

        domain_x = 2.0_real64 * pi
        domain_y = 2.0_real64 * pi

        ! Periodic grid (excludes duplicate endpoint)
        dx2 = domain_x / nx
        dy2 = domain_y / ny

        allocate(x2d(nx), y2d(ny), f2d(nx,ny), df2d_exact(nx,ny), df2d_compact(nx,ny))

        ! Generate grid points
        do i = 1, nx
            x2d(i) = (i - 1) * dx2
        end do

        do j = 1, ny
            y2d(j) = (j - 1) * dy2
        end do

        ! Initialize function f(x,y) = sin(x)*cos(y)
        do j = 1, ny
            do i = 1, nx
                f2d(i,j) = sin(x2d(i)) * cos(y2d(j))
            end do
        end do

        ! Compute exact x-derivative: df/dx = cos(x)*cos(y)
        do j = 1, ny
            do i = 1, nx
                df2d_exact(i,j) = cos(x2d(i)) * cos(y2d(j))
            end do
        end do

        ! Compute x-derivative using 4th-order compact scheme (periodic)
        df2d_compact = compute_compact_derivative(f2d, dx2, dy2, 1, 1, lhs_4th_order, rhs_4th_order, periodic=.true.)

        ! Compute error metrics (all points, since periodic has no boundary issues)
        l2_error = sqrt(sum((df2d_compact - df2d_exact)**2) / (nx*ny))
        max_error = maxval(abs(df2d_compact - df2d_exact))

        ! Print results
        print *, "  Grid size: ", nx, "x", ny
        print *, "  Grid spacing: dx=", dx2, ", dy=", dy2
        print *, "  Compact scheme results for x-derivative:"
        print *, "    L2 error: ", l2_error
        print *, "    Max error: ", max_error
        print *, ""

        deallocate(x2d, y2d, f2d, df2d_exact, df2d_compact)
    end subroutine test_2d_derivative

    ! Subroutine to write results to a file for plotting
    subroutine write_results(filename, x, f, df_exact, df_compact, df_fd)
        character(len=*), intent(in) :: filename
        real(real64), intent(in) :: x(:), f(:), df_exact(:), df_compact(:), df_fd(:)
        integer :: i, unit

        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A)') "# x    f(x)    df_exact    df_compact    df_fd    error_compact    error_fd"

        do i = 1, size(x)
            write(unit, '(7ES16.8)') x(i), f(i), df_exact(i), df_compact(i), df_fd(i), &
                                    abs(df_compact(i) - df_exact(i)), abs(df_fd(i) - df_exact(i))
        end do

        close(unit)
    end subroutine write_results

end program test_compact_derivatives