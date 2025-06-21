module findDuplicateSymmetryOperators_m
    use constants, only: eps12

    implicit none

    public :: findDuplicateSymmetryOperators

    private

contains

subroutine findDuplicateSymmetryOperators( &
    lattice_vectors, symm_op, n_symm_op, debug &
)
    double precision, intent(in) :: lattice_vectors(3,3)
    integer, intent(inout) :: symm_op(3,4,n_symm_op), n_symm_op
    logical, intent(in) :: debug

    integer :: ix, jx, kx, iix, io, jo, counter
    double precision :: lattice_translation(3), delta_t(3)
    logical :: filter(n_symm_op)

    filter(:) = .true.

    do io = 1, n_symm_op
        if(.not. filter(io)) cycle
        do jo = 1, n_symm_op
            if(io .eq. jo) cycle
            if( &
                (.not. (symm_op(1,1,io) .eq. symm_op(1,1,jo))) .or. &
                (.not. (symm_op(1,2,io) .eq. symm_op(1,2,jo))) .or. &
                (.not. (symm_op(1,3,io) .eq. symm_op(1,3,jo))) .or. &
                (.not. (symm_op(2,1,io) .eq. symm_op(2,1,jo))) .or. &
                (.not. (symm_op(2,2,io) .eq. symm_op(2,2,jo))) .or. &
                (.not. (symm_op(2,3,io) .eq. symm_op(2,3,jo))) .or. &
                (.not. (symm_op(3,1,io) .eq. symm_op(3,1,jo))) .or. &
                (.not. (symm_op(3,2,io) .eq. symm_op(3,2,jo))) .or. &
                (.not. (symm_op(3,3,io) .eq. symm_op(3,3,jo))) &
            ) cycle

            do iix = 1, 3
                delta_t(iix) = symm_op(iix,4,io) - symm_op(iix,4,jo)
            enddo

            do ix = -12, 12
                do jx = -12, 12
                    do kx = -12, 12
                        do iix = 1, 3
                            lattice_translation(iix) = &
                                ix * lattice_vectors(1,iix) + &
                                jx * lattice_vectors(2,iix) + &
                                kx * lattice_vectors(3,iix)
                        enddo

                        if( &
                            (dble(delta_t(1))/12.d0 .eq. lattice_translation(1)) .and. &
                            (dble(delta_t(2))/12.d0 .eq. lattice_translation(2)) .and. &
                            (dble(delta_t(3))/12.d0 .eq. lattice_translation(3)) &
                        ) then
                            if(debug) then
                                write(6,'(a)') "findDuplicateSymmetryOperators: Equivalent translation found for"
                                write(6,'(a,i4,a,i4)') "findDuplicateSymmetryOperators: ", io, " and ", jo
                            endif

                            if( &
                                (symm_op(1,4,jo) .eq. 0) .and. &
                                (symm_op(2,4,jo) .eq. 0) .and. &
                                (symm_op(3,4,jo) .eq. 0) &
                            ) then
                                filter(io) = .false.
                            else
                                filter(jo) = .false.
                            endif
                        endif
                    enddo
                    if(.not. filter(jo)) exit
                enddo
                if(.not. filter(jo)) exit
            enddo
        enddo
    enddo

    counter = 0
    do io = 1, n_symm_op
        if(filter(io)) then
            counter = counter + 1
            symm_op(1:3,1:4,counter) = symm_op(1:3,1:4,io)
        endif
    enddo
    n_symm_op = count(filter(1:n_symm_op))

    if(debug) then
        write(6,'(a)') "findDuplicateSymmetryOperators: Symmetry operators after removing duplicates"
        write(6,'(a,i5)') "findDuplicateSymmetryOperators: Number of operators =", n_symm_op
        do io = 1, n_symm_op
            write(6,'(a,i5)') "findDuplicateSymmetryOperators: Index = ", io
            write(6,'(4i4)') symm_op(1,1:4,io)
            write(6,'(4i4)') symm_op(2,1:4,io)
            write(6,'(4i4)') symm_op(3,1:4,io)
        enddo
    endif

end subroutine findDuplicateSymmetryOperators

end module findDuplicateSymmetryOperators_m