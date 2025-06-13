module restricAtomicCoordinates_m
    use constants, only: eps16, eps6
    use calculatereciprocallattice_m, only: calculateReciprocalLattice

    implicit none

    public :: restricAtomicCoordinates

    private

contains

subroutine restricAtomicCoordinates( &
    lattice_vectors, n_atom, atomic_coordinates, debug &
)
    integer, intent(in) :: n_atom
    double precision, intent(in) :: lattice_vectors(3,3)
    double precision, intent(inout) :: atomic_coordinates(3,n_atom)
    logical, intent(in) :: debug

    integer :: ia, i
    double precision :: temp_atomic_coordinate(3), reciprocal_lattice_vectors(3,3)

    ! Change atomic coordinate formats to fractional 
    ! and make the interval -0.5 < x <= 0.5
    call calculateReciprocalLattice(lattice_vectors, reciprocal_lattice_vectors, 0)

    if(debug) write(6,*) "restricAtomicCoordinates: Atomic coordinates in fractional coordinate"
    if(debug) write(6,*) "restricAtomicCoordinates: with restriction -0.5 < x <= 0.5"
    do ia = 1, n_atom
        temp_atomic_coordinate(1:3) = atomic_coordinates(1:3,ia)
        do i = 1, 3
            atomic_coordinates(i,ia) = &
                reciprocal_lattice_vectors(i,1) * temp_atomic_coordinate(1) + &
                reciprocal_lattice_vectors(i,2) * temp_atomic_coordinate(2) + &
                reciprocal_lattice_vectors(i,3) * temp_atomic_coordinate(3)
            ! Makes the interval -0.5 < x <= 0.5
            atomic_coordinates(i,ia) = &
                modulo(atomic_coordinates(i,ia)+0.5d0,1.d0) - 0.5d0
            ! This is to take care of exactly 0.5d0 
            if( atomic_coordinates(i,ia) + 0.5d0 .le. eps16) then
                atomic_coordinates(i,ia) = 0.5d0
            endif
        enddo
        if(debug) write(6,'(3f16.9)') atomic_coordinates(1:3,ia)
    enddo

end subroutine restricAtomicCoordinates

end module restricAtomicCoordinates_m