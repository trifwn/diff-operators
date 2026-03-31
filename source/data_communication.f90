module data_communication
    use iso_fortran_env, only: real64
    use MPI
    implicit none

    private
    public :: commit_array, free_vars, distribute, collect
    public :: exchange_halos

    ! MPI-related variables
    integer :: rank, size_mpi, ierr
    integer :: cart_comm

    ! DATA-related variables
    integer :: ndim
    integer, allocatable :: dims(:), coords(:)
    integer, allocatable :: global_sizes(:), local_sizes(:)
    integer, allocatable :: padded_counts(:), padded_starts(:, :), padded_sizes(:, :)
    integer, allocatable :: true_counts(:), true_starts(:, :), true_sizes(:, :)
    integer :: arr_padding = 1
    logical, allocatable :: periodic_dims(:)

    interface collect
        module subroutine collect_4D(global_data, local_data)
            real(real64), allocatable, intent(inout) :: global_data(:, :, :, :)
            real(real64), allocatable, intent(in) :: local_data(:, :, :, :)
        end subroutine collect_4D

        module subroutine collect_3D(global_data, local_data)
            real(real64), allocatable, intent(inout) :: global_data(:, :, :)
            real(real64), allocatable, intent(in) :: local_data(:, :, :)
        end subroutine collect_3D

        module subroutine collect_2D(global_data, local_data)
            real(real64), allocatable, intent(inout) :: global_data(:, :)
            real(real64), allocatable, intent(in) :: local_data(:, :)
        end subroutine collect_2D
    end interface collect

    interface distribute
        module subroutine distribute_4D(global_data, local_data)
            real(real64), target, intent(in) :: global_data(:, :, :, :)
            real(real64), allocatable, intent(out) :: local_data(:, :, :, :)
        end subroutine distribute_4D

        module subroutine distribute_3D(global_data, local_data)
            real(real64), target, intent(in) :: global_data(:, :, :)
            real(real64), allocatable, intent(out) :: local_data(:, :, :)
        end subroutine distribute_3D

        module subroutine distribute_2D(global_data, local_data)
            real(real64), target, intent(in) :: global_data(:, :)
            real(real64), allocatable, intent(out) :: local_data(:, :)
        end subroutine distribute_2D
    end interface distribute

    interface exchange_halos
        module procedure exchange_halos_3d
        module procedure exchange_halos_2d
    end interface exchange_halos

contains

    subroutine check_mpi(caller, mpi_ierr)
        character(len=*), intent(in) :: caller
        integer, intent(in) :: mpi_ierr
        if (mpi_ierr /= MPI_SUCCESS) then
            write(*, '(A,A,A,I0)') 'MPI error in ', caller, ': code=', mpi_ierr
            call MPI_Abort(cart_comm, mpi_ierr, ierr)
        end if
    end subroutine check_mpi

    subroutine free_vars()
        if (allocated(padded_counts)) deallocate(padded_counts)
        if (allocated(padded_starts)) deallocate(padded_starts)
        if (allocated(padded_sizes)) deallocate(padded_sizes)
        if (allocated(true_counts)) deallocate(true_counts)
        if (allocated(true_starts)) deallocate(true_starts)
        if (allocated(true_sizes)) deallocate(true_sizes)
        if (allocated(dims)) deallocate(dims)
        if (allocated(coords)) deallocate(coords)
        if (allocated(global_sizes)) deallocate(global_sizes)
        if (allocated(local_sizes)) deallocate(local_sizes)
        if (allocated(periodic_dims)) deallocate(periodic_dims)
        call MPI_Comm_free(cart_comm, ierr)
        ndim = 0
        arr_padding = 1
    end subroutine free_vars

    subroutine commit_array(array, dimensions, padding, periodic)
        implicit none
        real(real64), intent(in) :: array(..)
        integer, intent(in) :: dimensions(:)
        integer, intent(in) :: padding
        logical, intent(in), optional :: periodic(:)
        integer :: i

        ! Initialize MPI
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call check_mpi('commit_array:Comm_rank', ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size_mpi, ierr)
        call check_mpi('commit_array:Comm_size', ierr)

        arr_padding = padding

        ! Get dimensions of global data
        ndim = size(shape(array))
        if (allocated(global_sizes)) deallocate(global_sizes)
        if (allocated(local_sizes)) deallocate(local_sizes)
        if (allocated(periodic_dims)) deallocate(periodic_dims)
        allocate(global_sizes(ndim))
        allocate(local_sizes(ndim))
        allocate(periodic_dims(ndim))

        if (allocated(dims)) deallocate(dims)
        if (allocated(coords)) deallocate(coords)
        allocate(dims(ndim))
        allocate(coords(ndim))

        if (rank == 0) then
            do i = 1, ndim
                global_sizes(i) = size(array, i)
            end do
        end if
        ! Broadcast global sizes to all processes
        call MPI_Bcast(global_sizes, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call check_mpi('commit_array:Bcast', ierr)
        dims = dimensions

        ! Configurable periodicity
        if (present(periodic)) then
            periodic_dims = periodic
        else
            periodic_dims = .false.
        end if

        call MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periodic_dims, &
                             .true., cart_comm, ierr)
        call check_mpi('commit_array:Cart_create', ierr)
        call MPI_Cart_coords(cart_comm, rank, ndim, coords, ierr)
        call check_mpi('commit_array:Cart_coords', ierr)
        call calculate_access()
    end subroutine commit_array

    subroutine calculate_access()
        integer :: local_coords(ndim), local_size
        integer :: proc, i

        if (allocated(padded_counts)) deallocate(padded_counts)
        if (allocated(padded_starts)) deallocate(padded_starts)
        if (allocated(padded_sizes)) deallocate(padded_sizes)
        if (allocated(true_counts)) deallocate(true_counts)
        if (allocated(true_starts)) deallocate(true_starts)
        if (allocated(true_sizes)) deallocate(true_sizes)

        allocate(padded_counts(size_mpi))
        allocate(true_counts(size_mpi))

        allocate(padded_starts(size_mpi, ndim))
        allocate(true_starts(size_mpi, ndim))

        allocate(padded_sizes(size_mpi, ndim))
        allocate(true_sizes(size_mpi, ndim))

        padded_counts = 0
        padded_starts = 0
        do proc = 0, size_mpi - 1
            call MPI_Cart_coords(cart_comm, proc, ndim, local_coords, ierr)
            do i = 1, ndim
                local_size = global_sizes(i) / dims(i)
                if (local_coords(i) == dims(i) - 1) then
                    local_size = global_sizes(i) - local_size * (dims(i) - 1)
                    true_starts(proc + 1, i) = global_sizes(i) - local_size
                else
                    true_starts(proc + 1, i) = local_coords(i) * local_size
                end if
                true_sizes(proc + 1, i) = local_size

                ! Ghost cells: include padding from neighbors
                if (local_coords(i) > 0 .or. periodic_dims(i)) then
                    padded_starts(proc + 1, i) = true_starts(proc + 1, i) - arr_padding
                    ! Handle wrap-around for periodic
                    if (padded_starts(proc + 1, i) < 0) then
                        padded_starts(proc + 1, i) = padded_starts(proc + 1, i) + global_sizes(i)
                    end if
                    padded_sizes(proc + 1, i) = true_sizes(proc + 1, i) + arr_padding
                else
                    padded_starts(proc + 1, i) = true_starts(proc + 1, i)
                    padded_sizes(proc + 1, i) = true_sizes(proc + 1, i)
                end if

                if (local_coords(i) < dims(i) - 1 .or. periodic_dims(i)) then
                    padded_sizes(proc + 1, i) = padded_sizes(proc + 1, i) + arr_padding
                end if
            end do
            true_counts(proc + 1) = product(true_sizes(proc + 1, :))
            padded_counts(proc + 1) = product(padded_sizes(proc + 1, :))
            if (rank == proc) then
                local_sizes = true_sizes(proc + 1, :)
            end if
        end do
        call MPI_Barrier(cart_comm, ierr)
    end subroutine calculate_access

    ! -------------------------------------------------------------------------
    ! Distribute: root sends subdomains (with ghost cells) to all ranks
    ! -------------------------------------------------------------------------

    module subroutine distribute_4D(global_data, local_data)
        real(real64), target, intent(in) :: global_data(:, :, :, :)
        real(real64), allocatable, intent(out) :: local_data(:, :, :, :)
        integer :: source

        if (allocated(local_data)) deallocate(local_data)
        allocate(local_data(padded_sizes(rank + 1, 1), &
                            padded_sizes(rank + 1, 2), &
                            padded_sizes(rank + 1, 3), &
                            padded_sizes(rank + 1, 4)))
        source = 0
        if (rank == source) then
            call send_to_all(global_data)
            ! Root copies its own portion using padded_sizes (not local_sizes)
            local_data = global_data( &
                padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + padded_sizes(source + 1, 1), &
                padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + padded_sizes(source + 1, 2), &
                padded_starts(source + 1, 3) + 1:padded_starts(source + 1, 3) + padded_sizes(source + 1, 3), &
                padded_starts(source + 1, 4) + 1:padded_starts(source + 1, 4) + padded_sizes(source + 1, 4))
        else
            call recv_from_source(local_data, source)
        end if
        call MPI_Barrier(cart_comm, ierr)
    end subroutine distribute_4D

    module subroutine distribute_3D(global_data, local_data)
        real(real64), target, intent(in) :: global_data(:, :, :)
        real(real64), allocatable, intent(out) :: local_data(:, :, :)
        integer :: source

        if (allocated(local_data)) deallocate(local_data)
        allocate(local_data(padded_sizes(rank + 1, 1), &
                            padded_sizes(rank + 1, 2), &
                            padded_sizes(rank + 1, 3)))
        source = 0
        if (rank == source) then
            call send_to_all(global_data)
            local_data = global_data( &
                padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + padded_sizes(source + 1, 1), &
                padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + padded_sizes(source + 1, 2), &
                padded_starts(source + 1, 3) + 1:padded_starts(source + 1, 3) + padded_sizes(source + 1, 3))
        else
            call recv_from_source(local_data, source)
        end if
        call MPI_Barrier(cart_comm, ierr)
    end subroutine distribute_3D

    module subroutine distribute_2D(global_data, local_data)
        real(real64), target, intent(in) :: global_data(:, :)
        real(real64), allocatable, intent(out) :: local_data(:, :)
        integer :: source

        if (allocated(local_data)) deallocate(local_data)
        allocate(local_data(padded_sizes(rank + 1, 1), padded_sizes(rank + 1, 2)))
        source = 0
        if (rank == source) then
            call send_to_all(global_data)
            local_data = global_data( &
                padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + padded_sizes(source + 1, 1), &
                padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + padded_sizes(source + 1, 2))
        else
            call recv_from_source(local_data, source)
        end if
        call MPI_Barrier(cart_comm, ierr)
    end subroutine distribute_2D

    ! -------------------------------------------------------------------------
    ! Collect: all ranks send true (non-ghost) data back to root
    ! -------------------------------------------------------------------------

    module subroutine collect_4D(global_data, local_data)
        real(real64), allocatable, intent(inout) :: global_data(:, :, :, :)
        real(real64), allocatable, intent(in) :: local_data(:, :, :, :)
        integer :: source

        if (allocated(global_data)) deallocate(global_data)
        if (rank == 0) then
            allocate(global_data(global_sizes(1), global_sizes(2), &
                                 global_sizes(3), global_sizes(4)))
        end if

        source = 0
        if (rank == source) then
            call recv_from_all(global_data)
            global_data( &
                true_starts(source + 1, 1) + 1:true_starts(source + 1, 1) + true_sizes(source + 1, 1), &
                true_starts(source + 1, 2) + 1:true_starts(source + 1, 2) + true_sizes(source + 1, 2), &
                true_starts(source + 1, 3) + 1:true_starts(source + 1, 3) + true_sizes(source + 1, 3), &
                true_starts(source + 1, 4) + 1:true_starts(source + 1, 4) + true_sizes(source + 1, 4)) &
                = local_data( &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 1), &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 2), &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 3), &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 4))
        else
            call send_to_source(local_data, source)
        end if
    end subroutine collect_4D

    module subroutine collect_3D(global_data, local_data)
        real(real64), allocatable, intent(inout) :: global_data(:, :, :)
        real(real64), allocatable, intent(in) :: local_data(:, :, :)
        integer :: source

        if (allocated(global_data)) deallocate(global_data)
        if (rank == 0) then
            allocate(global_data(global_sizes(1), global_sizes(2), global_sizes(3)))
        end if

        source = 0
        if (rank == source) then
            call recv_from_all(global_data)
            global_data( &
                true_starts(source + 1, 1) + 1:true_starts(source + 1, 1) + true_sizes(source + 1, 1), &
                true_starts(source + 1, 2) + 1:true_starts(source + 1, 2) + true_sizes(source + 1, 2), &
                true_starts(source + 1, 3) + 1:true_starts(source + 1, 3) + true_sizes(source + 1, 3)) &
                = local_data( &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 1), &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 2), &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 3))
        else
            call send_to_source(local_data, source)
        end if
    end subroutine collect_3D

    module subroutine collect_2D(global_data, local_data)
        real(real64), allocatable, intent(inout) :: global_data(:, :)
        real(real64), allocatable, intent(in) :: local_data(:, :)
        integer :: source

        if (allocated(global_data)) deallocate(global_data)
        if (rank == 0) then
            allocate(global_data(global_sizes(1), global_sizes(2)))
        end if

        source = 0
        if (rank == source) then
            call recv_from_all(global_data)
            global_data( &
                true_starts(source + 1, 1) + 1:true_starts(source + 1, 1) + true_sizes(source + 1, 1), &
                true_starts(source + 1, 2) + 1:true_starts(source + 1, 2) + true_sizes(source + 1, 2)) &
                = local_data( &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 1), &
                    arr_padding + 1:arr_padding + true_sizes(source + 1, 2))
        else
            call send_to_source(local_data, source)
        end if
    end subroutine collect_2D

    ! -------------------------------------------------------------------------
    ! Halo exchange: exchange ghost cells between neighboring ranks
    ! -------------------------------------------------------------------------

    subroutine exchange_halos_3d(local_data)
        real(real64), intent(inout) :: local_data(:, :, :)
        integer :: nx, ny, nz
        integer :: src, dest, dim_idx
        integer :: status(MPI_STATUS_SIZE)
        real(real64), allocatable :: send_buf(:, :), recv_buf(:, :)

        nx = size(local_data, 1)
        ny = size(local_data, 2)
        nz = size(local_data, 3)

        ! Exchange along each dimension
        do dim_idx = 1, min(ndim, 3)
            call MPI_Cart_shift(cart_comm, dim_idx - 1, 1, src, dest, ierr)
            call check_mpi('exchange_halos_3d:Cart_shift', ierr)

            select case (dim_idx)
            case (1)
                ! Send right boundary, receive left ghost
                allocate(send_buf(ny, nz), recv_buf(ny, nz))
                ! Send to dest (right neighbor), recv from src (left neighbor)
                send_buf = local_data(nx - arr_padding, :, :)
                call MPI_Sendrecv(send_buf, ny * nz, MPI_DOUBLE_PRECISION, dest, 1, &
                                  recv_buf, ny * nz, MPI_DOUBLE_PRECISION, src, 1, &
                                  cart_comm, status, ierr)
                call check_mpi('exchange_halos_3d:Sendrecv_x_left', ierr)
                if (src /= MPI_PROC_NULL) local_data(1, :, :) = recv_buf

                ! Send left boundary, receive right ghost
                send_buf = local_data(1 + arr_padding, :, :)
                call MPI_Sendrecv(send_buf, ny * nz, MPI_DOUBLE_PRECISION, src, 2, &
                                  recv_buf, ny * nz, MPI_DOUBLE_PRECISION, dest, 2, &
                                  cart_comm, status, ierr)
                call check_mpi('exchange_halos_3d:Sendrecv_x_right', ierr)
                if (dest /= MPI_PROC_NULL) local_data(nx, :, :) = recv_buf
                deallocate(send_buf, recv_buf)

            case (2)
                allocate(send_buf(nx, nz), recv_buf(nx, nz))
                send_buf = local_data(:, ny - arr_padding, :)
                call MPI_Sendrecv(send_buf, nx * nz, MPI_DOUBLE_PRECISION, dest, 3, &
                                  recv_buf, nx * nz, MPI_DOUBLE_PRECISION, src, 3, &
                                  cart_comm, status, ierr)
                call check_mpi('exchange_halos_3d:Sendrecv_y_left', ierr)
                if (src /= MPI_PROC_NULL) local_data(:, 1, :) = recv_buf

                send_buf = local_data(:, 1 + arr_padding, :)
                call MPI_Sendrecv(send_buf, nx * nz, MPI_DOUBLE_PRECISION, src, 4, &
                                  recv_buf, nx * nz, MPI_DOUBLE_PRECISION, dest, 4, &
                                  cart_comm, status, ierr)
                call check_mpi('exchange_halos_3d:Sendrecv_y_right', ierr)
                if (dest /= MPI_PROC_NULL) local_data(:, ny, :) = recv_buf
                deallocate(send_buf, recv_buf)

            case (3)
                allocate(send_buf(nx, ny), recv_buf(nx, ny))
                send_buf = local_data(:, :, nz - arr_padding)
                call MPI_Sendrecv(send_buf, nx * ny, MPI_DOUBLE_PRECISION, dest, 5, &
                                  recv_buf, nx * ny, MPI_DOUBLE_PRECISION, src, 5, &
                                  cart_comm, status, ierr)
                call check_mpi('exchange_halos_3d:Sendrecv_z_left', ierr)
                if (src /= MPI_PROC_NULL) local_data(:, :, 1) = recv_buf

                send_buf = local_data(:, :, 1 + arr_padding)
                call MPI_Sendrecv(send_buf, nx * ny, MPI_DOUBLE_PRECISION, src, 6, &
                                  recv_buf, nx * ny, MPI_DOUBLE_PRECISION, dest, 6, &
                                  cart_comm, status, ierr)
                call check_mpi('exchange_halos_3d:Sendrecv_z_right', ierr)
                if (dest /= MPI_PROC_NULL) local_data(:, :, nz) = recv_buf
                deallocate(send_buf, recv_buf)
            end select
        end do
    end subroutine exchange_halos_3d

    subroutine exchange_halos_2d(local_data)
        real(real64), intent(inout) :: local_data(:, :)
        integer :: nx, ny
        integer :: src, dest, dim_idx
        integer :: status(MPI_STATUS_SIZE)
        real(real64), allocatable :: send_buf(:), recv_buf(:)

        nx = size(local_data, 1)
        ny = size(local_data, 2)

        do dim_idx = 1, min(ndim, 2)
            call MPI_Cart_shift(cart_comm, dim_idx - 1, 1, src, dest, ierr)
            call check_mpi('exchange_halos_2d:Cart_shift', ierr)

            select case (dim_idx)
            case (1)
                allocate(send_buf(ny), recv_buf(ny))
                send_buf = local_data(nx - arr_padding, :)
                call MPI_Sendrecv(send_buf, ny, MPI_DOUBLE_PRECISION, dest, 1, &
                                  recv_buf, ny, MPI_DOUBLE_PRECISION, src, 1, &
                                  cart_comm, status, ierr)
                if (src /= MPI_PROC_NULL) local_data(1, :) = recv_buf

                send_buf = local_data(1 + arr_padding, :)
                call MPI_Sendrecv(send_buf, ny, MPI_DOUBLE_PRECISION, src, 2, &
                                  recv_buf, ny, MPI_DOUBLE_PRECISION, dest, 2, &
                                  cart_comm, status, ierr)
                if (dest /= MPI_PROC_NULL) local_data(nx, :) = recv_buf
                deallocate(send_buf, recv_buf)

            case (2)
                allocate(send_buf(nx), recv_buf(nx))
                send_buf = local_data(:, ny - arr_padding)
                call MPI_Sendrecv(send_buf, nx, MPI_DOUBLE_PRECISION, dest, 3, &
                                  recv_buf, nx, MPI_DOUBLE_PRECISION, src, 3, &
                                  cart_comm, status, ierr)
                if (src /= MPI_PROC_NULL) local_data(:, 1) = recv_buf

                send_buf = local_data(:, 1 + arr_padding)
                call MPI_Sendrecv(send_buf, nx, MPI_DOUBLE_PRECISION, src, 4, &
                                  recv_buf, nx, MPI_DOUBLE_PRECISION, dest, 4, &
                                  cart_comm, status, ierr)
                if (dest /= MPI_PROC_NULL) local_data(:, ny) = recv_buf
                deallocate(send_buf, recv_buf)
            end select
        end do
    end subroutine exchange_halos_2d

    ! -------------------------------------------------------------------------
    ! Internal: point-to-point communication helpers
    ! Uses cart_comm consistently and non-blocking sends
    ! -------------------------------------------------------------------------

    subroutine recv_from_source(local_data, source)
        real(real64), intent(out) :: local_data(..)
        integer :: source, my_starts(ndim), my_size(ndim)
        integer :: status(MPI_STATUS_SIZE)
        integer :: recv_subarray

        my_starts = 0
        my_size = padded_sizes(rank + 1, :)
        call MPI_Type_create_subarray(ndim, my_size, my_size, my_starts, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      recv_subarray, ierr)
        call check_mpi('recv_from_source:Type_create', ierr)
        call MPI_Type_commit(recv_subarray, ierr)
        call MPI_Recv(local_data, 1, recv_subarray, source, 0, &
                      cart_comm, status, ierr)
        call check_mpi('recv_from_source:Recv', ierr)
        call MPI_Type_free(recv_subarray, ierr)
    end subroutine recv_from_source

    subroutine send_to_all(global_data)
        real(real64), intent(in) :: global_data(..)
        integer :: dest_starts(ndim), dest_size(ndim)
        integer :: i, num_sends
        integer, allocatable :: send_subarrays(:)
        integer, allocatable :: requests(:)
        integer, allocatable :: statuses(:, :)

        ! Count non-self ranks
        num_sends = size_mpi - 1
        if (num_sends == 0) return

        allocate(send_subarrays(num_sends))
        allocate(requests(num_sends))
        allocate(statuses(MPI_STATUS_SIZE, num_sends))

        num_sends = 0
        do i = 0, size_mpi - 1
            if (i == rank) cycle
            num_sends = num_sends + 1
            dest_starts = padded_starts(i + 1, :)
            dest_size = padded_sizes(i + 1, :)
            call MPI_Type_create_subarray(ndim, global_sizes, dest_size, dest_starts, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                          send_subarrays(num_sends), ierr)
            call check_mpi('send_to_all:Type_create', ierr)
            call MPI_Type_commit(send_subarrays(num_sends), ierr)
            call MPI_Isend(global_data, 1, send_subarrays(num_sends), i, 0, &
                           cart_comm, requests(num_sends), ierr)
            call check_mpi('send_to_all:Isend', ierr)
        end do

        call MPI_Waitall(num_sends, requests, statuses, ierr)
        call check_mpi('send_to_all:Waitall', ierr)

        ! Free types after all sends complete
        do i = 1, num_sends
            call MPI_Type_free(send_subarrays(i), ierr)
        end do

        deallocate(send_subarrays, requests, statuses)
    end subroutine send_to_all

    subroutine send_to_source(local_data, dest)
        real(real64), intent(in) :: local_data(..)
        integer, intent(in) :: dest
        integer :: send_subarray
        integer :: my_starts(ndim), my_sizes(ndim), my_true_sizes(ndim)

        my_starts = true_starts(rank + 1, :) - padded_starts(rank + 1, :)
        my_sizes = padded_sizes(rank + 1, :)
        my_true_sizes = true_sizes(rank + 1, :)
        call MPI_Type_create_subarray(ndim, my_sizes, my_true_sizes, my_starts, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      send_subarray, ierr)
        call check_mpi('send_to_source:Type_create', ierr)
        call MPI_Type_commit(send_subarray, ierr)
        call MPI_Send(local_data, 1, send_subarray, dest, 0, cart_comm, ierr)
        call check_mpi('send_to_source:Send', ierr)
        call MPI_Type_free(send_subarray, ierr)
    end subroutine send_to_source

    subroutine recv_from_all(global_data)
        real(real64), intent(out) :: global_data(..)
        integer :: dest, true_proc_starts(ndim), true_proc_sizes(ndim)
        integer :: num_recvs, i
        integer, allocatable :: recv_subarrays(:)
        integer, allocatable :: requests(:)
        integer, allocatable :: statuses(:, :)

        num_recvs = size_mpi - 1
        if (num_recvs == 0) return

        allocate(recv_subarrays(num_recvs))
        allocate(requests(num_recvs))
        allocate(statuses(MPI_STATUS_SIZE, num_recvs))

        num_recvs = 0
        do dest = 0, size_mpi - 1
            if (dest == rank) cycle
            num_recvs = num_recvs + 1
            true_proc_sizes = true_sizes(dest + 1, :)
            true_proc_starts = true_starts(dest + 1, :)
            call MPI_Type_create_subarray(ndim, global_sizes, true_proc_sizes, &
                                          true_proc_starts, MPI_ORDER_FORTRAN, &
                                          MPI_DOUBLE_PRECISION, &
                                          recv_subarrays(num_recvs), ierr)
            call check_mpi('recv_from_all:Type_create', ierr)
            call MPI_Type_commit(recv_subarrays(num_recvs), ierr)
            call MPI_Irecv(global_data, 1, recv_subarrays(num_recvs), dest, 0, &
                           cart_comm, requests(num_recvs), ierr)
            call check_mpi('recv_from_all:Irecv', ierr)
        end do

        call MPI_Waitall(num_recvs, requests, statuses, ierr)
        call check_mpi('recv_from_all:Waitall', ierr)

        do i = 1, num_recvs
            call MPI_Type_free(recv_subarrays(i), ierr)
        end do

        deallocate(recv_subarrays, requests, statuses)
    end subroutine recv_from_all

end module data_communication
