module symmetrize_m

    implicit none

    public :: symmetrize_vector, symmetrize_matrix

    private

contains

subroutine symmetrize_vector(n_symm_op, symm_op, tol, vector, debug)
    integer, intent(in) :: n_symm_op, symm_op(3,4,192)
    double precision, intent(in) :: tol
    double precision, intent(inout) :: vector(3)
    logical, intent(in) :: debug

    integer :: io, ix
    double precision :: &
        projected_vector(3), symmetric_vector(3)

    symmetric_vector(1:3) = 0.d0
    do io = 1, n_symm_op
        if( &
            (symm_op(1,1,io) .eq. 1) .and. &
            (symm_op(2,2,io) .eq. 1) .and. &
            (symm_op(3,3,io) .eq. 1) .and. &
            (symm_op(1,2,io) .eq. 0) .and. &
            (symm_op(2,1,io) .eq. 0) .and. &
            (symm_op(2,3,io) .eq. 0) .and. &
            (symm_op(3,2,io) .eq. 0) .and. &
            (symm_op(3,1,io) .eq. 0) .and. &
            (symm_op(1,3,io) .eq. 0) &
        ) cycle
        do ix = 1, 3
            projected_vector(ix) = &
                dble(symm_op(ix,1,io)) * vector(1) + &
                dble(symm_op(ix,2,io)) * vector(2) + &
                dble(symm_op(ix,3,io)) * vector(3) + &
                dble(symm_op(ix,4,io)) / 12.d0
            projected_vector(ix) = modulo(projected_vector(ix)+0.5d0,1.d0) - 0.5d0
        enddo

        if( &
            (abs(projected_vector(1) - vector(1)) .lt. tol) .and. &
            (abs(projected_vector(2) - vector(2)) .lt. tol) .and. &
            (abs(projected_vector(3) - vector(3)) .lt. tol) &
        ) &
            symmetric_vector(1:3) = &
                symmetric_vector(1:3) + projected_vector(1:3)
    enddo
    symmetric_vector(1:3) = symmetric_vector(1:3) / n_symm_op

    vector(1:3) = symmetric_vector(1:3)

end subroutine symmetrize_vector

subroutine symmetrize_matrix(n_symm_op, symm_op, matrix, debug)
    integer, intent(in) :: n_symm_op, symm_op(3,4,192)
    double precision, intent(inout) :: matrix(3,3)
    logical, intent(in) :: debug

    integer :: io, ix, jx
    double precision :: &
        temp_matrix(3,3), projected_matrix(3,3), symmetric_matrix(3,3)

    symmetric_matrix(1:3,1:3) = 0.d0
    do io = 1, n_symm_op
        temp_matrix(1:3,1:3) = matrix(1:3,1:3)
        do ix = 1, 3;do jx = 1, 3
            projected_matrix(ix,jx) = &
                symm_op(1,ix,io) * temp_matrix(jx,1) + &
                symm_op(2,ix,io) * temp_matrix(jx,2) + &
                symm_op(3,ix,io) * temp_matrix(jx,3)
        enddo;enddo
        if(debug) then
            write(6,'(a,i4)') "writeOutput: Projected matrix to operator index", io
            write(6,'(3f16.9)') projected_matrix(1,1:3)
            write(6,'(3f16.9)') projected_matrix(2,1:3)
            write(6,'(3f16.9)') projected_matrix(3,1:3)
        endif

        do ix = 1, 3;do jx = 1, 3
            symmetric_matrix(ix,jx) = &
                symmetric_matrix(ix,jx) + projected_matrix(ix,jx)
        enddo;enddo
    enddo
    symmetric_matrix(1:3,1:3) = &
        symmetric_matrix(1:3,1:3) / dble(n_symm_op)

    matrix(1:3,1:3) = symmetric_matrix(1:3,1:3)

end subroutine symmetrize_matrix

end module symmetrize_m