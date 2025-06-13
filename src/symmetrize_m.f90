module symmetrize_m

    implicit none

    public :: symmetrize_vector, symmetrize_matrix

    private

contains

subroutine symmetrize_vector(n_symm_op, symm_op, vector, debug)
    integer, intent(in) :: n_symm_op, symm_op(3,4,192)
    double precision, intent(inout) :: vector(3)
    logical, intent(in) :: debug

    integer :: io, ix
    double precision :: &
        projected_vector(3), symmetric_vector(3)

    symmetric_vector(1:3) = 0.d0
    do io = 1, n_symm_op
        do ix = 1, 3
            projected_vector(ix) = &
                symm_op(1,ix,io) * vector(1) + &
                symm_op(2,ix,io) * vector(2) + &
                symm_op(3,ix,io) * vector(3) - &
                symm_op(ix,4,io)
        enddo

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
                symm_op(1,ix,io) * temp_matrix(1,jx) + &
                symm_op(2,ix,io) * temp_matrix(2,jx) + &
                symm_op(3,ix,io) * temp_matrix(3,jx)
        enddo;enddo
        temp_matrix(1:3,1:3) = projected_matrix(1:3,1:3)
        do ix = 1, 3;do jx = 1, 3
            projected_matrix(ix,jx) = &
                temp_matrix(ix,1) * symm_op(1,jx,io) + &
                temp_matrix(ix,2) * symm_op(2,jx,io) + &
                temp_matrix(ix,3) * symm_op(3,jx,io)
        enddo;enddo

        symmetric_matrix(1:3,1:3) = &
            symmetric_matrix(1:3,1:3) + projected_matrix(1:3,1:3) 
    enddo
    symmetric_matrix(1:3,1:3) = &
        symmetric_matrix(1:3,1:3) / dble(n_symm_op)

    matrix(1:3,1:3) = symmetric_matrix(1:3,1:3)

end subroutine symmetrize_matrix

end module symmetrize_m