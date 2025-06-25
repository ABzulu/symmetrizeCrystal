module identifySpaceGroup_m
    implicit none

    public :: identifySpaceGroup

    private

contains

subroutine identifySpaceGroup( &
    symmetry_operators, n_symmetry_operators, &
    n_ops, rotations, translations, hall_symbol, found_space_group, debug &
)
    integer, intent(inout) :: &
        symmetry_operators(3,4,n_symmetry_operators), n_symmetry_operators, &
        n_ops(530), rotations(3,3,192,530), translations(3,192,530)
    character(len=17), intent(in) :: hall_symbol(530)
    logical, intent(in) :: debug
    logical, intent(out) :: found_space_group

    integer :: &
        ig, io, ix, jo, &
        n_shifts, n_denominators, denominators(11), &
        id, iid, jjd, kkd, shift(3), identity(3,3)
    integer, allocatable :: copy_symmetry_operators(:,:,:)
    logical :: w_found_in_group(192)

    found_space_group = .false.
    do ig = 1, 530
        w_found_in_group(:) = .false.

        if(.not. (n_ops(ig) .eq. n_symmetry_operators)) cycle

        do io = 1, n_symmetry_operators
            do jo = 1, n_ops(ig)
                if( &
                    (symmetry_operators(1,1,io) .eq. rotations(1,1,jo,ig))  .and. &
                    (symmetry_operators(1,2,io) .eq. rotations(1,2,jo,ig))  .and. &
                    (symmetry_operators(1,3,io) .eq. rotations(1,3,jo,ig))  .and. &
                    (symmetry_operators(2,1,io) .eq. rotations(2,1,jo,ig))  .and. &
                    (symmetry_operators(2,2,io) .eq. rotations(2,2,jo,ig))  .and. &
                    (symmetry_operators(2,3,io) .eq. rotations(2,3,jo,ig))  .and. &
                    (symmetry_operators(3,1,io) .eq. rotations(3,1,jo,ig))  .and. &
                    (symmetry_operators(3,2,io) .eq. rotations(3,2,jo,ig))  .and. &
                    (symmetry_operators(3,3,io) .eq. rotations(3,3,jo,ig))  .and. &
                    (symmetry_operators(1,4,io) .eq. translations(1,jo,ig)) .and. &
                    (symmetry_operators(2,4,io) .eq. translations(2,jo,ig)) .and. &
                    (symmetry_operators(3,4,io) .eq. translations(3,jo,ig)) &
                ) w_found_in_group(io) = .true.
            enddo

            if(w_found_in_group(io)) then
                cycle
            else
                exit
            endif
        enddo

        if(all(w_found_in_group(1:n_symmetry_operators))) then
            write(6,'(a)') "identifySpaceGroup: Successfully found Space group"
            write(6,'(a,a)') "identifySpaceGroup: Hall symbol = ", trim(hall_symbol(ig))
            found_space_group = .true.
            exit
        endif
    enddo

    if(debug .and. .not. found_space_group) then
        write(6,'(a)') "Could not find space group, trying shifts."
    else
        return
    endif

    n_denominators = 0
    denominators(:) = 0
    do io = 1, n_symmetry_operators
        do ix = 1, 3
            if(symmetry_operators(ix,4,io) .ne. 0) then
                if(n_denominators .eq. 0) then
                    n_denominators = n_denominators + 1
                    denominators(1) = symmetry_operators(ix,4,io)
                endif
                if(all(denominators(1:n_denominators) .ne. symmetry_operators(ix,4,io))) then
                    n_denominators = n_denominators + 1
                    denominators(n_denominators) = symmetry_operators(ix,4,io)
                endif
            endif
        enddo
    enddo

    n_shifts = 0
    do id = 1, n_denominators
        n_shifts = n_shifts + 12/denominators(id)*12/denominators(id)*12/denominators(id) - 1
    enddo

    if(debug) write(6,'(a,i8)') "identifySpaceGroup: Number of shifts to try = ", n_shifts

    allocate(copy_symmetry_operators(3,4,n_symmetry_operators))

    identity(1,1:3) = [ 1, 0, 0 ]
    identity(2,1:3) = [ 0, 1, 0 ]
    identity(3,1:3) = [ 0, 0, 1 ]
    copy_symmetry_operators(:,:,:) = symmetry_operators(:,:,:)
    do id = 1, n_denominators
        if(debug) &
            write(6,'(a,3i4)') "identifySpaceGroup: Denominator = ", 12/denominators(id)
        do iid = 0, 12/denominators(id) - 1;do jjd = 0, 12/denominators(id) - 1;do kkd = 0, 12/denominators(id) - 1 
            if((iid .eq. 0) .and. (jjd .eq. 0) .and. (kkd .eq. 0)) cycle
            shift(1) = iid * denominators(id)
            shift(2) = jjd * denominators(id)
            shift(3) = kkd * denominators(id)
            if(debug) &
                write(6,'(a,3i4)') "identifySpaceGroup: Trying origin shift = ", shift(1:3)
            do io = 1, n_symmetry_operators
                do ix = 1, 3
                    symmetry_operators(ix,4,io) = &
                        copy_symmetry_operators(ix,4,io) + &
                        copy_symmetry_operators(ix,1,io)*shift(1) + &
                        copy_symmetry_operators(ix,2,io)*shift(2) + &
                        copy_symmetry_operators(ix,3,io)*shift(3) - &
                        identity(ix,1)*shift(1) - &
                        identity(ix,2)*shift(2) - &
                        identity(ix,3)*shift(3)
                    if(symmetry_operators(ix,4,io) .ge. 12) &
                        symmetry_operators(ix,4,io) = &
                            symmetry_operators(ix,4,io) - 12
                    if(symmetry_operators(ix,4,io) .lt. 0) &
                        symmetry_operators(ix,4,io) = &
                            symmetry_operators(ix,4,io) + 12
                enddo

                if(debug) then
                    write(6,'(a,i4)') "identifySpaceGroup: Shifted operator, ", io
                    write(6,'(4i4)') symmetry_operators(1,1:4,io)
                    write(6,'(4i4)') symmetry_operators(2,1:4,io)
                    write(6,'(4i4)') symmetry_operators(3,1:4,io)
                endif
            enddo

            found_space_group = .false.
            do ig = 1, 530
                w_found_in_group(:) = .false.

                if(.not. (n_ops(ig) .eq. n_symmetry_operators)) cycle

                do io = 1, n_symmetry_operators
                    do jo = 1, n_ops(ig)
                        if( &
                            (symmetry_operators(1,1,io) .eq. rotations(1,1,jo,ig))  .and. &
                            (symmetry_operators(1,2,io) .eq. rotations(1,2,jo,ig))  .and. &
                            (symmetry_operators(1,3,io) .eq. rotations(1,3,jo,ig))  .and. &
                            (symmetry_operators(2,1,io) .eq. rotations(2,1,jo,ig))  .and. &
                            (symmetry_operators(2,2,io) .eq. rotations(2,2,jo,ig))  .and. &
                            (symmetry_operators(2,3,io) .eq. rotations(2,3,jo,ig))  .and. &
                            (symmetry_operators(3,1,io) .eq. rotations(3,1,jo,ig))  .and. &
                            (symmetry_operators(3,2,io) .eq. rotations(3,2,jo,ig))  .and. &
                            (symmetry_operators(3,3,io) .eq. rotations(3,3,jo,ig))  .and. &
                            (symmetry_operators(1,4,io) .eq. translations(1,jo,ig)) .and. &
                            (symmetry_operators(2,4,io) .eq. translations(2,jo,ig)) .and. &
                            (symmetry_operators(3,4,io) .eq. translations(3,jo,ig)) &
                        ) w_found_in_group(io) = .true.
                    enddo

                    if(w_found_in_group(io)) then
                        cycle
                    else
                        exit
                    endif
                enddo

                if(all(w_found_in_group(1:n_symmetry_operators))) then
                    write(6,'(a)') "identifySpaceGroup: Successfully found Space group"
                    write(6,'(a,a)') "identifySpaceGroup: Hall symbol = ", trim(hall_symbol(ig))
                    found_space_group = .true.
                    exit
                endif
            enddo
            if(found_space_group) then
                return
            else
                write(6,'(a)') "identifySpaceGroup: Shift did not result in space group"
            endif
        enddo;enddo;enddo
    enddo

    deallocate(copy_symmetry_operators)

end subroutine identifySpaceGroup

end module identifySpaceGroup_m