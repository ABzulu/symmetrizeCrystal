module identifySpaceGroup_m
    implicit none

    public :: identifySpaceGroup

    private

contains

subroutine identifySpaceGroup( &
    symm_op, n_symm_op, &
    n_ops, rotations, translations, hall_symbol, found_space_group, debug &
)
    integer, intent(inout) :: &
        symm_op(3,4,192), n_symm_op, &
        n_ops(530), rotations(3,3,192,530), translations(3,192,530)
    character(len=17), intent(in) :: hall_symbol(530)
    logical, intent(in) :: debug
    logical, intent(out) :: found_space_group

    integer :: ig, io, jo
    logical :: w_found_in_group(192)

    found_space_group = .false.
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
            found_space_group = .true.
            exit
        endif
    enddo

end subroutine identifySpaceGroup

end module identifySpaceGroup_m