program test_compact_derivatives
    use iso_fortran_env, only: real64, output_unit
    use compact_derivatives
    use custom_stencil_derivatives, only: compute_derivative
    implicit none

    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    integer :: n_passed, n_failed, n_total

    ! Compact scheme stencils (used across tests)
    ! 4th-order Pade first derivative (alpha=1/4)
    real(real64) :: lhs_4th(3), rhs_4th(5)
    ! 6th-order first derivative (alpha=1/3)
    real(real64) :: lhs_6th(3), rhs_6th(5)
    ! 4th-order second derivative (alpha=1/10, a=6/5)
    real(real64) :: lhs_d2(3), rhs_d2(5)

    ! Initialize stencils
    lhs_4th = [0.25_real64, 1.0_real64, 0.25_real64]
    rhs_4th = [0.0_real64, -0.75_real64, 0.0_real64, &
               0.75_real64, 0.0_real64]

    lhs_6th = [1.0_real64/3, 1.0_real64, 1.0_real64/3]
    rhs_6th = [-1.0_real64/36, -7.0_real64/9, 0.0_real64, &
               7.0_real64/9, 1.0_real64/36]

    lhs_d2 = [0.1_real64, 1.0_real64, 0.1_real64]
    rhs_d2 = [0.0_real64, 1.2_real64, -2.4_real64, &
              1.2_real64, 0.0_real64]

    n_passed = 0
    n_failed = 0
    n_total = 0

    print *, "================================================"
    print *, " Compact Derivatives Test Suite"
    print *, "================================================"
    print *, ""

    ! --- 1D Periodic Tests ---
    call test_1d_periodic_first_deriv()
    call test_1d_periodic_second_deriv()
    call test_1d_periodic_6th_order()
    call test_1d_periodic_high_wavenumber()

    ! --- 1D Non-Periodic Tests ---
    call test_1d_nonperiodic_first_deriv()
    call test_1d_nonperiodic_second_deriv()

    ! --- Convergence Rate Tests ---
    call test_convergence_4th_order()
    call test_convergence_6th_order()
    call test_convergence_2nd_deriv()

    ! --- 2D Tests ---
    call test_2d_x_derivative()
    call test_2d_y_derivative()
    call test_2d_second_derivative()

    ! --- 3D Tests ---
    call test_3d_derivatives()

    ! --- Vector Field Tests ---
    call test_vector_1d()
    call test_vector_2d()

    ! --- Summary ---
    print *, "================================================"
    print *, " RESULTS"
    print *, "================================================"
    print '(A,I3,A,I3,A)', "  ", n_passed, " / ", n_total, " tests PASSED"
    if (n_failed > 0) then
        print '(A,I3,A)', "  ", n_failed, " tests FAILED"
        stop 1
    else
        print *, "  All tests passed!"
    end if

contains

    ! ------------------------------------------------------------------
    ! Helper: check if error is below tolerance, print PASS/FAIL
    ! ------------------------------------------------------------------
    subroutine check(test_name, error, tolerance)
        character(len=*), intent(in) :: test_name
        real(real64), intent(in) :: error, tolerance

        n_total = n_total + 1
        if (error <= tolerance) then
            n_passed = n_passed + 1
            print '(A,A,A,ES10.3,A,ES10.3)', &
                "  PASS: ", test_name, &
                " (err=", error, ", tol=", tolerance
        else
            n_failed = n_failed + 1
            print '(A,A,A,ES10.3,A,ES10.3)', &
                "  FAIL: ", test_name, &
                " (err=", error, ", tol=", tolerance
        end if
    end subroutine check

    ! ------------------------------------------------------------------
    ! Helper: check convergence order
    ! ------------------------------------------------------------------
    subroutine check_order(test_name, order, min_order)
        character(len=*), intent(in) :: test_name
        real(real64), intent(in) :: order, min_order

        n_total = n_total + 1
        if (order >= min_order) then
            n_passed = n_passed + 1
            print '(A,A,A,F5.2,A,F5.2)', &
                "  PASS: ", test_name, &
                " (order=", order, ", min=", min_order
        else
            n_failed = n_failed + 1
            print '(A,A,A,F5.2,A,F5.2)', &
                "  FAIL: ", test_name, &
                " (order=", order, ", min=", min_order
        end if
    end subroutine check_order

    ! ------------------------------------------------------------------
    ! 1D Periodic: sin(x) first derivative
    ! ------------------------------------------------------------------
    subroutine test_1d_periodic_first_deriv()
        integer, parameter :: n = 100
        real(real64) :: dx, x(n), f(n), df_exact(n)
        real(real64), allocatable :: df(:)
        real(real64) :: max_err
        integer :: i

        print *, "--- 1D Periodic First Derivative ---"

        dx = 2.0_real64 * pi / n
        do i = 1, n
            x(i) = (i - 1) * dx
        end do
        f = sin(x)
        df_exact = cos(x)

        df = compute_compact_derivative(f, dx, 1, &
            lhs_4th, rhs_4th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("sin(x) d/dx, 4th order", max_err, 1.0e-7_real64)

        df = compute_compact_derivative(f, dx, 1, &
            lhs_6th, rhs_6th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("sin(x) d/dx, 6th order", max_err, 1.0e-9_real64)

        print *, ""
    end subroutine test_1d_periodic_first_deriv

    ! ------------------------------------------------------------------
    ! 1D Periodic: sin(x) second derivative
    ! ------------------------------------------------------------------
    subroutine test_1d_periodic_second_deriv()
        integer, parameter :: n = 100
        real(real64) :: dx, x(n), f(n), df_exact(n)
        real(real64), allocatable :: df(:)
        real(real64) :: max_err
        integer :: i

        print *, "--- 1D Periodic Second Derivative ---"

        dx = 2.0_real64 * pi / n
        do i = 1, n
            x(i) = (i - 1) * dx
        end do
        f = sin(x)
        df_exact = -sin(x)

        df = compute_compact_derivative(f, dx, 2, &
            lhs_d2, rhs_d2, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("sin(x) d2/dx2, 4th order", max_err, 1.0e-7_real64)

        print *, ""
    end subroutine test_1d_periodic_second_deriv

    ! ------------------------------------------------------------------
    ! 1D Periodic: exp(sin(x)) with 6th-order scheme
    ! ------------------------------------------------------------------
    subroutine test_1d_periodic_6th_order()
        integer, parameter :: n = 100
        real(real64) :: dx, x(n), f(n), df_exact(n)
        real(real64), allocatable :: df(:)
        real(real64) :: max_err
        integer :: i

        print *, "--- 1D Periodic 6th Order (complex function) ---"

        dx = 2.0_real64 * pi / n
        do i = 1, n
            x(i) = (i - 1) * dx
        end do
        f = exp(sin(x))
        df_exact = cos(x) * exp(sin(x))

        df = compute_compact_derivative(f, dx, 1, &
            lhs_6th, rhs_6th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("exp(sin(x)) d/dx, 6th order", max_err, 1.0e-8_real64)

        print *, ""
    end subroutine test_1d_periodic_6th_order

    ! ------------------------------------------------------------------
    ! 1D Periodic: high wavenumber sin(kx) resolution
    ! ------------------------------------------------------------------
    subroutine test_1d_periodic_high_wavenumber()
        integer, parameter :: n = 64
        real(real64) :: dx, x(n), f(n), df_exact(n)
        real(real64), allocatable :: df_compact(:), df_fd(:)
        real(real64) :: err_compact, err_fd
        integer :: i, k

        print *, "--- 1D Periodic High Wavenumber ---"

        dx = 2.0_real64 * pi / n
        do i = 1, n
            x(i) = (i - 1) * dx
        end do

        ! Test sin(kx) for k = 4, 8 -- compact should resolve
        ! much better than standard FD at moderate wavenumbers
        do k = 4, 8, 4
            f = sin(real(k, real64) * x)
            df_exact = real(k, real64) * cos(real(k, real64) * x)

            df_compact = compute_compact_derivative(f, dx, 1, &
                lhs_4th, rhs_4th, periodic=.true.)
            df_fd = compute_derivative(f, dx, 1, 1)

            err_compact = maxval(abs(df_compact - df_exact))
            err_fd = maxval(abs(df_fd - df_exact))

            call check("sin(kx) compact < FD, k=" // &
                char(ichar('0') + k), err_compact, err_fd)
        end do

        print *, ""
    end subroutine test_1d_periodic_high_wavenumber

    ! ------------------------------------------------------------------
    ! 1D Non-Periodic: x^3 first derivative
    ! ------------------------------------------------------------------
    subroutine test_1d_nonperiodic_first_deriv()
        integer, parameter :: n = 101
        real(real64) :: dx, x(n), f(n), df_exact(n)
        real(real64), allocatable :: df(:)
        real(real64) :: max_err_interior, max_err_boundary
        integer :: i

        print *, "--- 1D Non-Periodic First Derivative ---"

        dx = 2.0_real64 * pi / (n - 1)
        do i = 1, n
            x(i) = (i - 1) * dx
        end do

        ! x^3: Lele's 3rd-order closure is exact for polynomials up to degree 3
        f = x**3
        df_exact = 3.0_real64 * x**2
        df = compute_compact_derivative(f, dx, 1, lhs_4th, rhs_4th)
        max_err_interior = maxval(abs(df(3:n-2) - df_exact(3:n-2)))
        call check("x^3 d/dx interior", max_err_interior, 1.0e-10_real64)

        ! Boundary points should also be accurate (3rd order closure)
        max_err_boundary = max(abs(df(1) - df_exact(1)), &
                               abs(df(n) - df_exact(n)))
        call check("x^3 d/dx boundary", max_err_boundary, 1.0e-10_real64)

        ! exp(x): non-periodic, non-polynomial
        f = exp(x)
        df_exact = exp(x)
        df = compute_compact_derivative(f, dx, 1, lhs_4th, rhs_4th)
        max_err_interior = maxval(abs(df(3:n-2) - df_exact(3:n-2)))
        call check("exp(x) d/dx non-periodic", &
            max_err_interior, 2.0e-3_real64)

        print *, ""
    end subroutine test_1d_nonperiodic_first_deriv

    ! ------------------------------------------------------------------
    ! 1D Non-Periodic: second derivative
    ! ------------------------------------------------------------------
    subroutine test_1d_nonperiodic_second_deriv()
        integer, parameter :: n = 101
        real(real64) :: dx, x(n), f(n), df_exact(n)
        real(real64), allocatable :: df(:)
        real(real64) :: max_err
        integer :: i

        print *, "--- 1D Non-Periodic Second Derivative ---"

        dx = 2.0_real64 * pi / (n - 1)
        do i = 1, n
            x(i) = (i - 1) * dx
        end do

        ! x^4: f''(x) = 12x^2
        f = x**4
        df_exact = 12.0_real64 * x**2
        df = compute_compact_derivative(f, dx, 2, lhs_d2, rhs_d2)
        max_err = maxval(abs(df(3:n-2) - df_exact(3:n-2)))
        call check("x^4 d2/dx2 interior", max_err, 1.0e-3_real64)

        ! exp(x): f''(x) = exp(x)
        f = exp(x)
        df_exact = exp(x)
        df = compute_compact_derivative(f, dx, 2, lhs_d2, rhs_d2)
        max_err = maxval(abs(df(3:n-2) - df_exact(3:n-2)))
        call check("exp(x) d2/dx2 non-periodic", &
            max_err, 2.0e-2_real64)

        print *, ""
    end subroutine test_1d_nonperiodic_second_deriv

    ! ------------------------------------------------------------------
    ! Convergence: 4th-order first derivative
    ! ------------------------------------------------------------------
    subroutine test_convergence_4th_order()
        integer, parameter :: n_grids = 4
        integer :: ns(n_grids), ig, i, n
        real(real64) :: dx, max_err(n_grids), order
        real(real64), allocatable :: x(:), f(:), df_exact(:), df(:)

        print *, "--- Convergence: 4th Order First Derivative ---"

        ns = [32, 64, 128, 256]

        do ig = 1, n_grids
            n = ns(ig)
            dx = 2.0_real64 * pi / n
            allocate(x(n), f(n), df_exact(n))
            do i = 1, n
                x(i) = (i - 1) * dx
            end do
            f = sin(x)
            df_exact = cos(x)
            df = compute_compact_derivative(f, dx, 1, &
                lhs_4th, rhs_4th, periodic=.true.)
            max_err(ig) = maxval(abs(df - df_exact))
            deallocate(x, f, df_exact)
        end do

        ! Compute convergence orders between successive grids
        do ig = 2, n_grids
            order = log(max_err(ig-1) / max_err(ig)) / log(2.0_real64)
            print '(A,I4,A,I4,A,F6.2)', &
                "    n=", ns(ig-1), " -> ", ns(ig), &
                ": order = ", order
        end do

        ! Check the last refinement gives ~4th order
        order = log(max_err(n_grids-1) / max_err(n_grids)) &
                / log(2.0_real64)
        call check_order("4th-order convergence rate", order, 3.8_real64)

        print *, ""
    end subroutine test_convergence_4th_order

    ! ------------------------------------------------------------------
    ! Convergence: 6th-order first derivative
    ! ------------------------------------------------------------------
    subroutine test_convergence_6th_order()
        integer, parameter :: n_grids = 4
        integer :: ns(n_grids), ig, i, n
        real(real64) :: dx, max_err(n_grids), order
        real(real64), allocatable :: x(:), f(:), df_exact(:), df(:)

        print *, "--- Convergence: 6th Order First Derivative ---"

        ns = [32, 64, 128, 256]

        do ig = 1, n_grids
            n = ns(ig)
            dx = 2.0_real64 * pi / n
            allocate(x(n), f(n), df_exact(n))
            do i = 1, n
                x(i) = (i - 1) * dx
            end do
            f = sin(x)
            df_exact = cos(x)
            df = compute_compact_derivative(f, dx, 1, &
                lhs_6th, rhs_6th, periodic=.true.)
            max_err(ig) = maxval(abs(df - df_exact))
            deallocate(x, f, df_exact)
        end do

        do ig = 2, n_grids
            order = log(max_err(ig-1) / max_err(ig)) / log(2.0_real64)
            print '(A,I4,A,I4,A,F6.2)', &
                "    n=", ns(ig-1), " -> ", ns(ig), &
                ": order = ", order
        end do

        order = log(max_err(n_grids-1) / max_err(n_grids)) &
                / log(2.0_real64)
        call check_order("6th-order convergence rate", order, 5.8_real64)

        print *, ""
    end subroutine test_convergence_6th_order

    ! ------------------------------------------------------------------
    ! Convergence: 4th-order second derivative
    ! ------------------------------------------------------------------
    subroutine test_convergence_2nd_deriv()
        integer, parameter :: n_grids = 4
        integer :: ns(n_grids), ig, i, n
        real(real64) :: dx, max_err(n_grids), order
        real(real64), allocatable :: x(:), f(:), df_exact(:), df(:)

        print *, "--- Convergence: 4th Order Second Derivative ---"

        ns = [32, 64, 128, 256]

        do ig = 1, n_grids
            n = ns(ig)
            dx = 2.0_real64 * pi / n
            allocate(x(n), f(n), df_exact(n))
            do i = 1, n
                x(i) = (i - 1) * dx
            end do
            f = sin(x)
            df_exact = -sin(x)
            df = compute_compact_derivative(f, dx, 2, &
                lhs_d2, rhs_d2, periodic=.true.)
            max_err(ig) = maxval(abs(df - df_exact))
            deallocate(x, f, df_exact)
        end do

        do ig = 2, n_grids
            order = log(max_err(ig-1) / max_err(ig)) / log(2.0_real64)
            print '(A,I4,A,I4,A,F6.2)', &
                "    n=", ns(ig-1), " -> ", ns(ig), &
                ": order = ", order
        end do

        order = log(max_err(n_grids-1) / max_err(n_grids)) &
                / log(2.0_real64)
        call check_order("d2/dx2 4th-order convergence", &
            order, 3.8_real64)

        print *, ""
    end subroutine test_convergence_2nd_deriv

    ! ------------------------------------------------------------------
    ! 2D: x-derivative (periodic)
    ! ------------------------------------------------------------------
    subroutine test_2d_x_derivative()
        integer, parameter :: nx = 50, ny = 50
        real(real64) :: dx, dy
        real(real64) :: x(nx), y(ny)
        real(real64) :: f(nx,ny), df_exact(nx,ny)
        real(real64), allocatable :: df(:,:)
        real(real64) :: max_err
        integer :: i, j

        print *, "--- 2D x-derivative (periodic) ---"

        dx = 2.0_real64 * pi / nx
        dy = 2.0_real64 * pi / ny
        do i = 1, nx
            x(i) = (i - 1) * dx
        end do
        do j = 1, ny
            y(j) = (j - 1) * dy
        end do

        do j = 1, ny
            do i = 1, nx
                f(i,j) = sin(x(i)) * cos(y(j))
                df_exact(i,j) = cos(x(i)) * cos(y(j))
            end do
        end do

        df = compute_compact_derivative(f, dx, dy, 1, 1, &
            lhs_4th, rhs_4th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("2D df/dx sin(x)*cos(y)", max_err, 1.0e-5_real64)

        print *, ""
    end subroutine test_2d_x_derivative

    ! ------------------------------------------------------------------
    ! 2D: y-derivative (periodic)
    ! ------------------------------------------------------------------
    subroutine test_2d_y_derivative()
        integer, parameter :: nx = 50, ny = 50
        real(real64) :: dx, dy
        real(real64) :: x(nx), y(ny)
        real(real64) :: f(nx,ny), df_exact(nx,ny)
        real(real64), allocatable :: df(:,:)
        real(real64) :: max_err
        integer :: i, j

        print *, "--- 2D y-derivative (periodic) ---"

        dx = 2.0_real64 * pi / nx
        dy = 2.0_real64 * pi / ny
        do i = 1, nx
            x(i) = (i - 1) * dx
        end do
        do j = 1, ny
            y(j) = (j - 1) * dy
        end do

        do j = 1, ny
            do i = 1, nx
                f(i,j) = sin(x(i)) * cos(y(j))
                df_exact(i,j) = -sin(x(i)) * sin(y(j))
            end do
        end do

        df = compute_compact_derivative(f, dx, dy, 2, 1, &
            lhs_4th, rhs_4th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("2D df/dy sin(x)*cos(y)", max_err, 1.0e-5_real64)

        print *, ""
    end subroutine test_2d_y_derivative

    ! ------------------------------------------------------------------
    ! 2D: second derivative (periodic)
    ! ------------------------------------------------------------------
    subroutine test_2d_second_derivative()
        integer, parameter :: nx = 50, ny = 50
        real(real64) :: dx, dy
        real(real64) :: x(nx), y(ny)
        real(real64) :: f(nx,ny), df_exact(nx,ny)
        real(real64), allocatable :: df(:,:)
        real(real64) :: max_err
        integer :: i, j

        print *, "--- 2D Second Derivative (periodic) ---"

        dx = 2.0_real64 * pi / nx
        dy = 2.0_real64 * pi / ny
        do i = 1, nx
            x(i) = (i - 1) * dx
        end do
        do j = 1, ny
            y(j) = (j - 1) * dy
        end do

        do j = 1, ny
            do i = 1, nx
                f(i,j) = sin(x(i)) * cos(y(j))
            end do
        end do

        ! d2f/dx2 = -sin(x)*cos(y)
        do j = 1, ny
            do i = 1, nx
                df_exact(i,j) = -sin(x(i)) * cos(y(j))
            end do
        end do

        df = compute_compact_derivative(f, dx, dy, 1, 2, &
            lhs_d2, rhs_d2, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("2D d2f/dx2 sin(x)*cos(y)", max_err, 1.0e-5_real64)

        ! d2f/dy2 = -sin(x)*cos(y)
        df = compute_compact_derivative(f, dx, dy, 2, 2, &
            lhs_d2, rhs_d2, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("2D d2f/dy2 sin(x)*cos(y)", max_err, 1.0e-5_real64)

        print *, ""
    end subroutine test_2d_second_derivative

    ! ------------------------------------------------------------------
    ! 3D: derivatives in all dimensions (periodic)
    ! ------------------------------------------------------------------
    subroutine test_3d_derivatives()
        integer, parameter :: nx = 20, ny = 20, nz = 20
        real(real64) :: dx, dy, dz
        real(real64) :: x(nx), y(ny), z(nz)
        real(real64) :: f(nx,ny,nz), df_exact(nx,ny,nz)
        real(real64), allocatable :: df(:,:,:)
        real(real64) :: max_err
        integer :: i, j, k

        print *, "--- 3D Derivatives (periodic) ---"

        dx = 2.0_real64 * pi / nx
        dy = 2.0_real64 * pi / ny
        dz = 2.0_real64 * pi / nz
        do i = 1, nx
            x(i) = (i - 1) * dx
        end do
        do j = 1, ny
            y(j) = (j - 1) * dy
        end do
        do k = 1, nz
            z(k) = (k - 1) * dz
        end do

        ! f(x,y,z) = sin(x) * cos(y) * sin(z)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    f(i,j,k) = sin(x(i)) * cos(y(j)) * sin(z(k))
                end do
            end do
        end do

        ! df/dx = cos(x) * cos(y) * sin(z)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    df_exact(i,j,k) = cos(x(i)) * cos(y(j)) * sin(z(k))
                end do
            end do
        end do
        df = compute_compact_derivative(f, dx, dy, dz, 1, 1, &
            lhs_4th, rhs_4th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("3D df/dx", max_err, 1.0e-3_real64)

        ! df/dy = -sin(x) * sin(y) * sin(z)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    df_exact(i,j,k) = -sin(x(i)) * sin(y(j)) * sin(z(k))
                end do
            end do
        end do
        df = compute_compact_derivative(f, dx, dy, dz, 2, 1, &
            lhs_4th, rhs_4th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("3D df/dy", max_err, 1.0e-3_real64)

        ! df/dz = sin(x) * cos(y) * cos(z)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    df_exact(i,j,k) = sin(x(i)) * cos(y(j)) * cos(z(k))
                end do
            end do
        end do
        df = compute_compact_derivative(f, dx, dy, dz, 3, 1, &
            lhs_4th, rhs_4th, periodic=.true.)
        max_err = maxval(abs(df - df_exact))
        call check("3D df/dz", max_err, 1.0e-3_real64)

        print *, ""
    end subroutine test_3d_derivatives

    ! ------------------------------------------------------------------
    ! 1D Vector field derivative (periodic)
    ! ------------------------------------------------------------------
    subroutine test_vector_1d()
        integer, parameter :: n = 64, neq = 3
        real(real64) :: dx
        real(real64) :: x(n)
        real(real64) :: f(neq, n), df_exact(neq, n)
        real(real64), allocatable :: df(:,:)
        real(real64) :: max_err
        integer :: i

        print *, "--- 1D Vector Field Derivative (periodic) ---"

        dx = 2.0_real64 * pi / n
        do i = 1, n
            x(i) = (i - 1) * dx
        end do

        ! 3-component vector: [sin(x), cos(x), sin(2x)]
        f(1,:) = sin(x)
        f(2,:) = cos(x)
        f(3,:) = sin(2.0_real64 * x)

        df_exact(1,:) = cos(x)
        df_exact(2,:) = -sin(x)
        df_exact(3,:) = 2.0_real64 * cos(2.0_real64 * x)

        df = compute_compact_derivative(f, dx, 1, 1, &
            lhs_4th, rhs_4th, periodic=.true.)

        max_err = maxval(abs(df(1,:) - df_exact(1,:)))
        call check("vec1D comp1: d(sin)/dx", max_err, 1.0e-6_real64)

        max_err = maxval(abs(df(2,:) - df_exact(2,:)))
        call check("vec1D comp2: d(cos)/dx", max_err, 1.0e-6_real64)

        max_err = maxval(abs(df(3,:) - df_exact(3,:)))
        call check("vec1D comp3: d(sin2x)/dx", max_err, 2.0e-5_real64)

        print *, ""
    end subroutine test_vector_1d

    ! ------------------------------------------------------------------
    ! 2D Vector field derivative (periodic)
    ! ------------------------------------------------------------------
    subroutine test_vector_2d()
        integer, parameter :: nx = 32, ny = 32, neq = 2
        real(real64) :: dx, dy
        real(real64) :: x(nx), y(ny)
        real(real64) :: f(neq, nx, ny), df_exact(neq, nx, ny)
        real(real64), allocatable :: df(:,:,:)
        real(real64) :: max_err
        integer :: i, j

        print *, "--- 2D Vector Field Derivative (periodic) ---"

        dx = 2.0_real64 * pi / nx
        dy = 2.0_real64 * pi / ny
        do i = 1, nx
            x(i) = (i - 1) * dx
        end do
        do j = 1, ny
            y(j) = (j - 1) * dy
        end do

        ! 2-component vector field:
        !   comp1: sin(x)*cos(y)
        !   comp2: cos(x)*sin(y)
        do j = 1, ny
            do i = 1, nx
                f(1, i, j) = sin(x(i)) * cos(y(j))
                f(2, i, j) = cos(x(i)) * sin(y(j))
            end do
        end do

        ! Test x-derivative (dim=1)
        ! d/dx comp1 = cos(x)*cos(y)
        ! d/dx comp2 = -sin(x)*sin(y)
        do j = 1, ny
            do i = 1, nx
                df_exact(1, i, j) = cos(x(i)) * cos(y(j))
                df_exact(2, i, j) = -sin(x(i)) * sin(y(j))
            end do
        end do

        df = compute_compact_derivative(f, dx, dy, 1, 1, &
            lhs_4th, rhs_4th, periodic=.true.)

        max_err = maxval(abs(df(1,:,:) - df_exact(1,:,:)))
        call check("vec2D comp1 df/dx", max_err, 1.0e-4_real64)

        max_err = maxval(abs(df(2,:,:) - df_exact(2,:,:)))
        call check("vec2D comp2 df/dx", max_err, 1.0e-4_real64)

        ! Test y-derivative (dim=2)
        ! d/dy comp1 = -sin(x)*sin(y)
        ! d/dy comp2 = cos(x)*cos(y)
        do j = 1, ny
            do i = 1, nx
                df_exact(1, i, j) = -sin(x(i)) * sin(y(j))
                df_exact(2, i, j) = cos(x(i)) * cos(y(j))
            end do
        end do

        df = compute_compact_derivative(f, dx, dy, 2, 1, &
            lhs_4th, rhs_4th, periodic=.true.)

        max_err = maxval(abs(df(1,:,:) - df_exact(1,:,:)))
        call check("vec2D comp1 df/dy", max_err, 1.0e-4_real64)

        max_err = maxval(abs(df(2,:,:) - df_exact(2,:,:)))
        call check("vec2D comp2 df/dy", max_err, 1.0e-4_real64)

        print *, ""
    end subroutine test_vector_2d

end program test_compact_derivatives