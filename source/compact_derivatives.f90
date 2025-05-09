module compact_derivatives
    use iso_fortran_env, only: real64
    implicit none
    
    private
    public :: compute_compact_derivative 
    
    ! Interfaces for the main compact derivative computation functions
    interface compute_compact_derivative
        module function compute_compact_derivative_1d(f, dx, order, lhs_stencil, rhs_stencil) result(df)
            implicit none
            real(real64), intent(in) :: f(:)
            real(real64), intent(in) :: dx
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)  ! [alpha_{i-1}, alpha_i, alpha_{i+1}]
            real(real64), intent(in) :: rhs_stencil(:)  ! [a_i-2, a_i-1, a_i, a_i+1, a_i+2]
            real(real64), allocatable :: df(:)
        end function
        
        module function compute_compact_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:)
            real(real64), intent(in) :: dx, dy
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            real(real64), allocatable :: df(:,:)
        end function
        
        module function compute_compact_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:,:)
            real(real64), intent(in) :: dx, dy, dz
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            real(real64), allocatable :: df(:,:,:)
        end function
        
        module function compute_compact_vector_derivative_1d(f, dx, dim, order, lhs_stencil, rhs_stencil) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:)
            real(real64), intent(in) :: dx
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            real(real64), allocatable :: df(:,:)
        end function
        
        module function compute_compact_vector_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:,:)
            real(real64), intent(in) :: dx, dy
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            real(real64), allocatable :: df(:,:,:)
        end function
        
        module function compute_compact_vector_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil) result(df)
            implicit none
            real(real64), intent(in) :: f(:,:,:,:)
            real(real64), intent(in) :: dx, dy, dz
            integer, intent(in) :: dim
            integer, intent(in) :: order
            real(real64), intent(in) :: lhs_stencil(:)
            real(real64), intent(in) :: rhs_stencil(:)
            real(real64), allocatable :: df(:,:,:,:)
        end function
    end interface compute_compact_derivative

    ! Private helper functions
    private :: tdma_solver

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
    
    ! Helper function to apply the RHS stencil for a specific point
    function apply_rhs_stencil(f, i, n, rhs_stencil, dx, order) result(rhs_val)
        real(real64), intent(in) :: f(:)
        integer, intent(in) :: i, n
        real(real64), intent(in) :: rhs_stencil(:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: order
        real(real64) :: rhs_val
        
        integer :: j, idx, stencil_size, offset
        
        stencil_size = size(rhs_stencil)
        offset = stencil_size / 2
        
        rhs_val = 0.0_real64
        
        do j = 1, stencil_size
            idx = i + (j - offset - 1)
            
            ! Apply boundary conditions for points outside the domain
            if (idx < 1) then
                idx = 1  ! Simple extrapolation to boundary
            else if (idx > n) then
                idx = n  ! Simple extrapolation to boundary
            end if
            
            rhs_val = rhs_val + rhs_stencil(j) * f(idx)
        end do
        
        ! Scale by appropriate power of dx
        rhs_val = rhs_val / (dx**order)
    end function apply_rhs_stencil
    
    ! 1D compact derivative computation
    module function compute_compact_derivative_1d(f, dx, order, lhs_stencil, rhs_stencil) result(df)
        real(real64), intent(in) :: f(:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)  ! [alpha_{i-1}, alpha_i, alpha_{i+1}]
        real(real64), intent(in) :: rhs_stencil(:)  ! [a_i-2, a_i-1, a_i, a_i+1, a_i+2]
        real(real64), allocatable :: df(:)
        
        integer :: n, i
        real(real64), allocatable :: a(:), b(:), c(:), d(:)
        
        n = size(f)
        allocate(df(n))
        allocate(a(n), b(n), c(n), d(n))
        
        ! Set up tridiagonal system coefficients
        ! For most interior points
        do i = 2, n-1
            ! LHS coefficients
            a(i) = lhs_stencil(1)  ! coefficient for f'_{i-1}
            b(i) = lhs_stencil(2)  ! coefficient for f'_i
            c(i) = lhs_stencil(3)  ! coefficient for f'_{i+1}
            
            ! RHS computation
            d(i) = apply_rhs_stencil(f, i, n, rhs_stencil, dx, order)
        end do
        
        ! Special treatment for boundaries (assuming 3-point LHS stencil)
        ! Left boundary (i=1): use one-sided approximation
        a(1) = 0.0_real64                ! No f'_{0} term
        b(1) = lhs_stencil(2)            ! Coefficient for f'_1
        c(1) = lhs_stencil(3)            ! Coefficient for f'_2
        d(1) = apply_rhs_stencil(f, 1, n, rhs_stencil, dx, order)
        
        ! Right boundary (i=n): use one-sided approximation
        a(n) = lhs_stencil(1)            ! Coefficient for f'_{n-1}
        b(n) = lhs_stencil(2)            ! Coefficient for f'_n
        c(n) = 0.0_real64                ! No f'_{n+1} term
        d(n) = apply_rhs_stencil(f, n, n, rhs_stencil, dx, order)
        
        ! Solve the tridiagonal system
        call tdma_solver(a, b, c, d, df, n)
    end function compute_compact_derivative_1d
    
    ! 2D compact derivative computation
    module function compute_compact_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil) result(df)
        real(real64), intent(in) :: f(:,:)
        real(real64), intent(in) :: dx, dy
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
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
                temp = compute_compact_derivative_1d(f(:,j), dx, order, lhs_stencil, rhs_stencil)
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
                temp = compute_compact_derivative_1d(f(i,:), dy, order, lhs_stencil, rhs_stencil)
                df(i,:) = temp
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel
        end select
    end function compute_compact_derivative_2d
    
    ! 3D compact derivative computation
    module function compute_compact_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil) result(df)
        real(real64), intent(in) :: f(:,:,:)
        real(real64), intent(in) :: dx, dy, dz
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
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
                    temp = compute_compact_derivative_1d(f(:,j,k), dx, order, lhs_stencil, rhs_stencil)
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
                    temp = compute_compact_derivative_1d(f(i,:,k), dy, order, lhs_stencil, rhs_stencil)
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
                    temp = compute_compact_derivative_1d(f(i,j,:), dz, order, lhs_stencil, rhs_stencil)
                    df(i,j,:) = temp
                end do
            end do
            !$omp end do
            deallocate(temp)
            !$omp end parallel
        end select
    end function compute_compact_derivative_3d
    
    ! Vector version of 1D compact derivative
    module function compute_compact_vector_derivative_1d(f, dx, dim, order, lhs_stencil, rhs_stencil) result(df)
        real(real64), intent(in) :: f(:,:)
        real(real64), intent(in) :: dx
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        real(real64), allocatable :: df(:,:)
        
        integer :: neq, n, i
        
        neq = size(f, 1)
        n = size(f, 2)
        allocate(df(neq, n))
        
        !$omp parallel do
        do i = 1, neq
            df(i,:) = compute_compact_derivative_1d(f(i,:), dx, order, lhs_stencil, rhs_stencil)
        end do
        !$omp end parallel do
    end function compute_compact_vector_derivative_1d
    
    ! Vector version of 2D compact derivative
    module function compute_compact_vector_derivative_2d(f, dx, dy, dim, order, lhs_stencil, rhs_stencil) result(df)
        real(real64), intent(in) :: f(:,:,:)
        real(real64), intent(in) :: dx, dy
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        real(real64), allocatable :: df(:,:,:)
        
        integer :: neq, nx, ny, i
        
        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        allocate(df(neq, nx, ny))
        
        !$omp parallel do
        do i = 1, neq
            df(i,:,:) = compute_compact_derivative_2d(f(i,:,:), dx, dy, dim, order, lhs_stencil, rhs_stencil)
        end do
        !$omp end parallel do
    end function compute_compact_vector_derivative_2d
    
    ! Vector version of 3D compact derivative
    module function compute_compact_vector_derivative_3d(f, dx, dy, dz, dim, order, lhs_stencil, rhs_stencil) result(df)
        real(real64), intent(in) :: f(:,:,:,:)
        real(real64), intent(in) :: dx, dy, dz
        integer, intent(in) :: dim
        integer, intent(in) :: order
        real(real64), intent(in) :: lhs_stencil(:)
        real(real64), intent(in) :: rhs_stencil(:)
        real(real64), allocatable :: df(:,:,:,:)
        
        integer :: neq, nx, ny, nz, i
        
        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)
        allocate(df(neq, nx, ny, nz))
        
        !$omp parallel do
        do i = 1, neq
            df(i,:,:,:) = compute_compact_derivative_3d(f(i,:,:,:), dx, dy, dz, dim, order, lhs_stencil, rhs_stencil)
        end do
        !$omp end parallel do
    end function compute_compact_vector_derivative_3d
    
end module compact_derivatives