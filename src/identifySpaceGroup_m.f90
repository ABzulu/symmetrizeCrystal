module identifySpaceGroup_m
    implicit none

    public :: identifySpaceGroup

    private

contains

subroutine identifySpaceGroup( &
    symm_op, n_symm_op, &
    n_ops, rotations, translations, hall_symbol, debug &
)
    integer, intent(inout) :: &
        symm_op(3,4,192), n_symm_op, &
        n_ops(530), rotations(3,3,192,530), translations(3,192,530)
    character(len=17), intent(in) :: hall_symbol(530)
    logical, intent(in) :: debug

    integer :: ig, io, jo, ix, jx, temp_symm_op(3,4), M(3,3)
    logical :: w_found_in_group(192)

    ! do io = 1, n_symm_op
    !     temp_symm_op(:,:) = symm_op(:,:,io)
    !     do ix = 1, 3; do jx = 1, 3
    !         symm_op(ix,jx,io) = &
    !             M(1,ix) * temp_symm_op(1,jx) + &
    !             M(2,ix) * temp_symm_op(2,jx) + &
    !             M(3,ix) * temp_symm_op(3,jx)
    !     enddo;enddo
    !     temp_symm_op(:,:) = symm_op(:,:,io)
    !     do ix = 1, 3; do jx = 1, 3
    !         symm_op(ix,jx,io) = &
    !             temp_symm_op(ix,1) * M(1,jx) + &
    !             temp_symm_op(ix,2) * M(2,jx) + &
    !             temp_symm_op(ix,3) * M(3,jx)
    !     enddo;enddo
    !     do ix = 1, 3
    !         symm_op(ix,4,io) = &
    !             M(1,ix) * temp_symm_op(ix,4) + &
    !             M(2,ix) * temp_symm_op(ix,4) + &
    !             M(3,ix) * temp_symm_op(ix,4)
    !     enddo
    ! enddo

    if(debug) then
        write(6,'(a)') "identifySpaceGroup: After changing basis to convenctional cell"
        do io = 1, n_symm_op
            write(6,'(a,i4)') "identifySpaceGroup: Symmetry operator", io
            do ix = 1, 3
                write(6,'(4i4)') symm_op(ix,1:4,io)
            enddo
        enddo
    endif

    do ig = 1, 530
        w_found_in_group(:) = .false.

        if(.not. (n_ops(ig) .eq. n_symm_op)) cycle

        do io = 1, n_symm_op
            do jo = 1, n_ops(ig)
                if( &
                    (symm_op(1,1,io) .eq. rotations(1,1,jo,ig))  .and. &
                    (symm_op(1,2,io) .eq. rotations(1,2,jo,ig))  .and. &
                    (symm_op(1,3,io) .eq. rotations(1,3,jo,ig))  .and. &
                    (symm_op(2,1,io) .eq. rotations(2,1,jo,ig))  .and. &
                    (symm_op(2,2,io) .eq. rotations(2,2,jo,ig))  .and. &
                    (symm_op(2,3,io) .eq. rotations(2,3,jo,ig))  .and. &
                    (symm_op(3,1,io) .eq. rotations(3,1,jo,ig))  .and. &
                    (symm_op(3,2,io) .eq. rotations(3,2,jo,ig))  .and. &
                    (symm_op(3,3,io) .eq. rotations(3,3,jo,ig))  .and. &
                    (symm_op(1,4,io) .eq. translations(1,jo,ig)) .and. &
                    (symm_op(2,4,io) .eq. translations(2,jo,ig)) .and. &
                    (symm_op(3,4,io) .eq. translations(3,jo,ig)) &
                ) w_found_in_group(io) = .true.
            enddo

            if(w_found_in_group(io)) then
                cycle
            else
                exit
            endif
        enddo

        if(all(w_found_in_group(1:n_symm_op))) then
            write(6,'(a)') "identifySpaceGroup: Successfully found Space group"
            write(6,'(a,a)') "identifySpaceGroup: Hall symbol = ", trim(hall_symbol(ig))
            exit
        endif
    enddo

    if(all(w_found_in_group(1:n_symm_op))) then
        continue
    else
        write(6,'(a)') "identifySpaceGroup: Could not find matching group"
        stop
    endif

end subroutine identifySpaceGroup

end module identifySpaceGroup_m