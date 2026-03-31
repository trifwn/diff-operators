program bench_scaling
    use iso_fortran_env, only: real64, int64
    use compact_derivatives
    use custom_stencil_derivatives, only: compute_derivative
    implicit none

    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)

    ! Stencils
    real(real64) :: lhs_4th(3), rhs_4th(5)
    real(real64) :: lhs_d2(3), rhs_d2(5)

    lhs_4th = [0.25_real64, 1.0_real64, 0.25_real64]
    rhs_4th = [0.0_real64, -0.75_real64, 0.0_real64, &
               0.75_real64, 0.0_real64]
    lhs_d2 = [0.1_real64, 1.0_real64, 0.1_real64]
    rhs_d2 = [0.0_real64, 1.2_real64, -2.4_real64, &
              1.2_real64, 0.0_real64]

    print *, "========================================================"
    print *, " Compact Derivatives: Scaling Benchmark"
    print *, "========================================================"
    print *, ""

    ! --- CPU Scaling: problem size ---
    call bench_1d_size_scaling()
    call bench_2d_size_scaling()
    call bench_3d_size_scaling()

    ! --- CPU Scaling: thread count ---
    call bench_thread_scaling()

    ! --- Memory analysis ---
    call bench_memory_scaling()

contains

    ! ------------------------------------------------------------------
    ! Timing helper
    ! ------------------------------------------------------------------
    subroutine get_time(t)
        real(real64), intent(out) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, real64) / real(rate, real64)
    end subroutine get_time

    ! ------------------------------------------------------------------
    ! 1D: time vs problem size
    ! ------------------------------------------------------------------
    subroutine bench_1d_size_scaling()
        integer, parameter :: n_sizes = 6, n_reps = 50
        integer :: sizes(n_sizes)
        integer :: is, ir, i, n
        real(real64) :: dx, t0, t1
        real(real64) :: time_compact, time_fd
        real(real64), allocatable :: x(:), f(:)
        real(real64), allocatable :: df(:)

        sizes = [100, 500, 1000, 5000, 10000, 50000]

        print *, "--- 1D Size Scaling ---"
        print '(A6, A14, A14, A10)', &
            "N", "Compact(us)", "FD(us)", "Ratio"

        do is = 1, n_sizes
            n = sizes(is)
            dx = 2.0_real64 * pi / n
            allocate(x(n), f(n))
            do i = 1, n
                x(i) = (i - 1) * dx
            end do
            f = sin(x) + 0.5_real64 * cos(3.0_real64 * x)

            ! Time compact scheme
            call get_time(t0)
            do ir = 1, n_reps
                df = compute_compact_derivative(f, dx, 1, &
                    lhs_4th, rhs_4th, periodic=.true.)
            end do
            call get_time(t1)
            time_compact = (t1 - t0) / n_reps * 1.0e6_real64

            ! Time standard FD
            call get_time(t0)
            do ir = 1, n_reps
                df = compute_derivative(f, dx, 1, 1)
            end do
            call get_time(t1)
            time_fd = (t1 - t0) / n_reps * 1.0e6_real64

            print '(I6, F14.1, F14.1, F10.2)', &
                n, time_compact, time_fd, time_compact / time_fd

            deallocate(x, f)
        end do
        print *, ""
    end subroutine bench_1d_size_scaling

    ! ------------------------------------------------------------------
    ! 2D: time vs problem size
    ! ------------------------------------------------------------------
    subroutine bench_2d_size_scaling()
        integer, parameter :: n_sizes = 5, n_reps = 10
        integer :: sizes(n_sizes)
        integer :: is, ir, i, j, n
        real(real64) :: dx, dy, t0, t1
        real(real64) :: time_compact, time_fd
        real(real64), allocatable :: f2d(:,:), df2d(:,:)
        real(real64), allocatable :: x(:), y(:)

        sizes = [50, 100, 200, 400, 800]

        print *, "--- 2D Size Scaling (N x N grid) ---"
        print '(A6, A10, A14, A14, A10)', &
            "N", "N^2", "Compact(us)", "FD(us)", "Ratio"

        do is = 1, n_sizes
            n = sizes(is)
            dx = 2.0_real64 * pi / n
            dy = dx
            allocate(x(n), y(n), f2d(n,n))
            do i = 1, n
                x(i) = (i - 1) * dx
                y(i) = (i - 1) * dy
            end do
            do j = 1, n
                do i = 1, n
                    f2d(i,j) = sin(x(i)) * cos(y(j))
                end do
            end do

            call get_time(t0)
            do ir = 1, n_reps
                df2d = compute_compact_derivative(f2d, dx, dy, &
                    1, 1, lhs_4th, rhs_4th, periodic=.true.)
            end do
            call get_time(t1)
            time_compact = (t1 - t0) / n_reps * 1.0e6_real64

            call get_time(t0)
            do ir = 1, n_reps
                df2d = compute_derivative(f2d, dx, dy, 1, 1)
            end do
            call get_time(t1)
            time_fd = (t1 - t0) / n_reps * 1.0e6_real64

            print '(I6, I10, F14.1, F14.1, F10.2)', &
                n, n*n, time_compact, time_fd, time_compact / time_fd

            deallocate(x, y, f2d)
        end do
        print *, ""
    end subroutine bench_2d_size_scaling

    ! ------------------------------------------------------------------
    ! 3D: time vs problem size
    ! ------------------------------------------------------------------
    subroutine bench_3d_size_scaling()
        integer, parameter :: n_sizes = 4, n_reps = 3
        integer :: sizes(n_sizes)
        integer :: is, ir, i, j, k, n
        real(real64) :: dx, dy, dz, t0, t1
        real(real64) :: time_compact, time_fd
        real(real64), allocatable :: f3d(:,:,:), df3d(:,:,:)
        real(real64), allocatable :: x(:), y(:), z(:)

        sizes = [20, 40, 80, 128]

        print *, "--- 3D Size Scaling (N x N x N grid) ---"
        print '(A6, A10, A14, A14, A10)', &
            "N", "N^3", "Compact(us)", "FD(us)", "Ratio"

        do is = 1, n_sizes
            n = sizes(is)
            dx = 2.0_real64 * pi / n
            dy = dx
            dz = dx
            allocate(x(n), y(n), z(n), f3d(n,n,n))
            do i = 1, n
                x(i) = (i - 1) * dx
                y(i) = (i - 1) * dy
                z(i) = (i - 1) * dz
            end do
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        f3d(i,j,k) = sin(x(i)) * cos(y(j)) * sin(z(k))
                    end do
                end do
            end do

            call get_time(t0)
            do ir = 1, n_reps
                df3d = compute_compact_derivative(f3d, dx, dy, dz, &
                    1, 1, lhs_4th, rhs_4th, periodic=.true.)
            end do
            call get_time(t1)
            time_compact = (t1 - t0) / n_reps * 1.0e6_real64

            call get_time(t0)
            do ir = 1, n_reps
                df3d = compute_derivative(f3d, dx, dy, dz, 1, 1)
            end do
            call get_time(t1)
            time_fd = (t1 - t0) / n_reps * 1.0e6_real64

            print '(I6, I10, F14.1, F14.1, F10.2)', &
                n, n*n*n, time_compact, time_fd, time_compact / time_fd

            deallocate(x, y, z, f3d)
        end do
        print *, ""
    end subroutine bench_3d_size_scaling

    ! ------------------------------------------------------------------
    ! Thread scaling: vary OMP_NUM_THREADS for a fixed 3D problem
    ! ------------------------------------------------------------------
    subroutine bench_thread_scaling()
        !$ use omp_lib
        integer, parameter :: n = 128, n_reps = 5
        integer :: ir, i, j, k, nt
        integer :: thread_counts(4)
        real(real64) :: dx, dy, dz, t0, t1
        real(real64) :: time_compact, time_base
        real(real64), allocatable :: f3d(:,:,:), df3d(:,:,:)
        real(real64) :: x(n), y(n), z(n)

        thread_counts = [1, 2, 3, 4]

        print *, "--- Thread Scaling (3D, N=128) ---"

        dx = 2.0_real64 * pi / n
        dy = dx
        dz = dx
        do i = 1, n
            x(i) = (i - 1) * dx
            y(i) = (i - 1) * dy
            z(i) = (i - 1) * dz
        end do
        allocate(f3d(n,n,n))
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    f3d(i,j,k) = sin(x(i)) * cos(y(j)) * sin(z(k))
                end do
            end do
        end do

        time_base = 0.0_real64

        print '(A8, A14, A12, A12)', &
            "Threads", "Time(us)", "Speedup", "Efficiency"

        do nt = 1, size(thread_counts)
            !$ call omp_set_num_threads(thread_counts(nt))

            ! Warm-up run
            df3d = compute_compact_derivative(f3d, dx, dy, dz, &
                1, 1, lhs_4th, rhs_4th, periodic=.true.)

            call get_time(t0)
            do ir = 1, n_reps
                df3d = compute_compact_derivative(f3d, dx, dy, dz, &
                    1, 1, lhs_4th, rhs_4th, periodic=.true.)
            end do
            call get_time(t1)
            time_compact = (t1 - t0) / n_reps * 1.0e6_real64

            if (nt == 1) time_base = time_compact

            print '(I8, F14.1, F12.2, F12.1, A1)', &
                thread_counts(nt), time_compact, &
                time_base / time_compact, &
                100.0_real64 * time_base / &
                    (time_compact * thread_counts(nt)), '%'
        end do

        deallocate(f3d)
        print *, ""
    end subroutine bench_thread_scaling

    ! ------------------------------------------------------------------
    ! Memory scaling analysis
    ! ------------------------------------------------------------------
    subroutine bench_memory_scaling()
        integer, parameter :: n_sizes = 5
        integer :: sizes(n_sizes)
        integer :: is, n
        real(real64) :: mem_input, mem_output, mem_tdma
        real(real64) :: mem_cyclic_extra, mem_total
        real(real64) :: bytes_per_real

        sizes = [64, 128, 256, 512, 1024]
        bytes_per_real = 8.0_real64  ! real64

        print *, "--- Memory Scaling Analysis ---"
        print *, ""

        ! 1D memory analysis
        print *, "1D Compact Derivative Memory (bytes):"
        print '(A8, A12, A12, A12, A12, A14)', &
            "N", "Input", "Output", "TDMA", "Cyclic", "Total"
        do is = 1, n_sizes
            n = sizes(is)
            mem_input = n * bytes_per_real                  ! f(n)
            mem_output = n * bytes_per_real                 ! df(n)
            mem_tdma = 4.0_real64 * n * bytes_per_real      ! a,b,c,d
            ! cyclic_tdma: bb(n), u(n), y(n), z(n) + cp(n), dp(n) from inner TDMA x2
            mem_cyclic_extra = 8.0_real64 * n * bytes_per_real
            mem_total = mem_input + mem_output + mem_tdma + mem_cyclic_extra
            print '(I8, F12.0, F12.0, F12.0, F12.0, F14.0)', &
                n, mem_input, mem_output, mem_tdma, &
                mem_cyclic_extra, mem_total
        end do
        print *, ""

        ! 3D memory analysis (per-thread for parallel regions)
        print *, "3D Compact Derivative Memory (MB, N x N x N grid):"
        print '(A8, A12, A12, A16, A16)', &
            "N", "Input(MB)", "Output(MB)", "PerThread(MB)", "Total4T(MB)"
        do is = 1, n_sizes
            n = sizes(is)
            if (n > 512) cycle  ! Skip very large 3D grids
            mem_input = real(n, real64)**3 * bytes_per_real / 1.0e6_real64
            mem_output = mem_input
            ! Per-thread: temp(n) array + 1D TDMA internals
            ! TDMA: a(n),b(n),c(n),d(n),cp(n),dp(n) = 6n
            ! Cyclic: bb(n),u(n),y(n),z(n) + 2*TDMA(cp,dp) = 4n + 4n = 8n
            ! temp(n) + df_1d(n) from function result = 2n
            ! Total per-thread: (6+8+2)*n = 16*n
            mem_tdma = 16.0_real64 * n * bytes_per_real / 1.0e6_real64
            ! Total with 4 threads
            mem_total = mem_input + mem_output + &
                        4.0_real64 * mem_tdma
            print '(I8, F12.2, F12.2, F16.4, F16.2)', &
                n, mem_input, mem_output, mem_tdma, mem_total
        end do
        print *, ""

        ! Summary
        print *, "--- Memory Scaling Summary ---"
        print *, "  1D: O(N) total memory"
        print *, "       - TDMA solver: 4N reals (a,b,c,d)"
        print *, "       - Cyclic TDMA: +8N reals (Sherman-Morrison)"
        print *, "       - Total: ~14N reals per solve"
        print *, ""
        print *, "  2D (NxM): O(N*M) for input/output arrays"
        print *, "       - Per-thread: O(max(N,M)) for 1D TDMA"
        print *, "       - Total parallel overhead: O(T*max(N,M))"
        print *, "       where T = number of threads"
        print *, ""
        print *, "  3D (NxNxN): O(N^3) for input/output arrays"
        print *, "       - Per-thread: O(N) for 1D TDMA"
        print *, "       - Total parallel overhead: O(T*N)"
        print *, "       - Overhead ratio: O(T/N^2) -> negligible"
        print *, ""
    end subroutine bench_memory_scaling

end program bench_scaling