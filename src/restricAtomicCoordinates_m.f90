module restricAtomicCoordinates_m
    use constants, only: eps16, eps6
    use calculatereciprocallattice_m, only: calculateReciprocalLattice

    implicit none

    public :: restricAtomicCoordinates

    private

contains

subroutine restricAtomicCoordinates( &
    modifyed_lattice_vectors, original_lattice_vectors, n_atom, atomic_coordinates, debug &
)
    integer, intent(in) :: n_atom
    double precision, intent(in) :: &
        modifyed_lattice_vectors(3,3), original_lattice_vectors(3,3)
    double precision, intent(inout) :: atomic_coordinates(3,n_atom)
    logical, intent(in) :: debug

    integer :: ia, i, ix, jx
    double precision :: &
        temp_atomic_coordinate(3), reciprocal_lattice_vectors(3,3), T(3,3)

    call calculateReciprocalLattice(original_lattice_vectors, reciprocal_lattice_vectors, 0)

    do ix = 1, 3;do jx = 1, 3
        T(ix,jx) = &
            reciprocal_lattice_vectors(ix,1) * modifyed_lattice_vectors(jx,1) + &
            reciprocal_lattice_vectors(ix,2) * modifyed_lattice_vectors(jx,2) + &
            reciprocal_lattice_vectors(ix,3) * modifyed_lattice_vectors(jx,3)
    enddo;enddo
    if(debug) then
        write(6,'(a)') "restrictAtomicCoordinates: T"
        write(6,'(3f16.9)') T(1,1:3)
        write(6,'(3f16.9)') T(2,1:3)
        write(6,'(3f16.9)') T(3,1:3)
    endif
    do ia = 1, n_atom
        temp_atomic_coordinate(1:3) = atomic_coordinates(1:3,ia)
        do ix = 1, 3
            atomic_coordinates(ix,ia) = &
                T(ix,1) * temp_atomic_coordinate(1) + &
                T(ix,2) * temp_atomic_coordinate(2) + &
                T(ix,3) * temp_atomic_coordinate(3)
        enddo
    enddo
    ! Change atomic coordinate formats to fractional 
    ! and make the interval -0.5 < x <= 0.5
    call calculateReciprocalLattice(modifyed_lattice_vectors, reciprocal_lattice_vectors, 0)

    if(debug) then
        write(6,'(a)') "restricAtomicCoordinates: Atomic coordinates before restriction"
        do ia = 1, n_atom
            write(6,'(3f16.9)') atomic_coordinates(1:3,ia)
        enddo
    endif

    if(debug) write(6,'(a)') "restricAtomicCoordinates: Atomic coordinates in fractional coordinate"
    if(debug) write(6,'(a)') "restricAtomicCoordinates: with restriction -0.5 < x <= 0.5"
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