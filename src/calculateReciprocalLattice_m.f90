module calculatereciprocallattice_m

    implicit none

    public :: calculateReciprocalLattice

    private

contains

subroutine calculateReciprocalLattice(A, B, iopt)
    ! Calculate reciprocal lattice vectors.
    ! Their product with direct lattice vectors is 1 if iopt = 0 or 2pi if iopt = 1
    ! AB^T = I not AB = I. This is to match the formating of A.

    implicit none

    double precision, intent(in)  :: A(3,3) 
    integer, intent(in) :: iopt
    double precision, intent(out) :: B(3,3)

    integer :: i
    double precision :: pi, C, CI

    pi = acos(-1.d0)

    B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
    B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
    B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
    B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
    B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
    B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
    B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
    B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)

    C = 1.d0
    if(iopt .eq. 1) C = 2.d0 * pi

    do i = 1,3
        CI=C/(A(1,i)*B(1,i)+A(2,i)*B(2,i)+A(3,i)*B(3,i))
        B(1,i)=B(1,i)*CI
        B(2,i)=B(2,i)*CI
        B(3,i)=B(3,i)*CI
    enddo

end subroutine calculateReciprocalLattice

end module calculatereciprocallattice_m