module structureinit_m
    use constants, only: eps16
    use calculatereciprocallattice_m, only: calculateReciprocalLattice

    implicit none

    public :: delaunayReduction, restricAtomicCoordinates

    private

contains

subroutine restricAtomicCoordinates( &
    lattice_vectors, n_atom, atomic_coordinates &
)
    integer, intent(in) :: n_atom
    double precision, intent(in) :: lattice_vectors(3,3)
    double precision, intent(inout) :: atomic_coordinates(3,n_atom)

    integer :: ia, i
    double precision :: temp_atomic_coordinate(3), reciprocal_lattice_vectors(3,3)

    ! Change atomic coordinate formats to fractional 
    ! and make the interval -0.5 < x <= 0.5
    call calculateReciprocalLattice(lattice_vectors, reciprocal_lattice_vectors, 0)
    do ia = 1, n_atom
        temp_atomic_coordinate(1:3) = atomic_coordinates(1:3,ia)
        do i = 1, 3
            atomic_coordinates(i,ia) = &
                reciprocal_lattice_vectors(1,i) * temp_atomic_coordinate(1) + &
                reciprocal_lattice_vectors(2,i) * temp_atomic_coordinate(2) + &
                reciprocal_lattice_vectors(3,i) * temp_atomic_coordinate(3)
            ! Makes the interval -0.5 < x <= 0.5
            atomic_coordinates(i,ia) = &
                modulo(atomic_coordinates(i,ia)+0.5d0,1.d0) - 0.5d0
            ! This is to take care of exactly 0.5d0 
            if( atomic_coordinates(i,ia) + 0.5d0 .le. eps16) then
                atomic_coordinates(i,ia) = 0.5d0
            endif
        enddo
    enddo

end subroutine restricAtomicCoordinates

subroutine delaunayReduction( &
    lattice_vectors, n_atom, atomic_coordinates, lattice_index &
)
    integer, intent(in) :: n_atom
    double precision, intent(inout) :: lattice_vectors(3,3), atomic_coordinates(3,n_atom)
    integer, intent(out) :: lattice_index(3)

    integer :: temp_i
    double precision :: a_sqr(3), temp_a, temp_lattice_vectors(3,3), scalar_product(3)
    double precision, allocatable :: temp_atomic_coordinates(:,:)
    logical :: done

    integer :: ir, i

    lattice_index(1) = 1
    lattice_index(2) = 2
    lattice_index(3) = 3
    do ir = 1, 100
        ! Order lattice vectors so that a_1 <= a_2 <= a_3
        do i = 1, 3
            a_sqr(i) = &
                lattice_vectors(1,lattice_index(i))*lattice_vectors(1,lattice_index(i)) + &
                lattice_vectors(2,lattice_index(i))*lattice_vectors(2,lattice_index(i)) + &
                lattice_vectors(3,lattice_index(i))*lattice_vectors(3,lattice_index(i))
        enddo
        if(a_sqr(2) .gt. a_sqr(3)) then
            temp_a = a_sqr(2)
            a_sqr(2) = a_sqr(3)
            a_sqr(3) = temp_a
            temp_i = lattice_index(2)
            lattice_index(2) = lattice_index(3)
            lattice_index(3) = temp_i
        endif
        if(a_sqr(1) .gt. a_sqr(2)) then
            temp_a = a_sqr(1)
            a_sqr(1) = a_sqr(2)
            a_sqr(2) = temp_a
            temp_i = lattice_index(1)
            lattice_index(1) = lattice_index(2)
            lattice_index(2) = temp_i
        endif
        if(a_sqr(2) .gt. a_sqr(3)) then
            temp_a = a_sqr(2)
            a_sqr(2) = a_sqr(3)
            a_sqr(3) = temp_a
            temp_i = lattice_index(2)
            lattice_index(2) = lattice_index(3)
            lattice_index(3) = temp_i
        endif

        ! Make sure scalar product of lattice vectors are <= 0
        scalar_product(1) = &
            lattice_vectors(1,lattice_index(1)) * lattice_vectors(1,lattice_index(2)) + &
            lattice_vectors(2,lattice_index(1)) * lattice_vectors(2,lattice_index(2)) + &
            lattice_vectors(3,lattice_index(1)) * lattice_vectors(3,lattice_index(2))
        scalar_product(2) = &
            lattice_vectors(1,lattice_index(2)) * lattice_vectors(1,lattice_index(3)) + &
            lattice_vectors(2,lattice_index(2)) * lattice_vectors(2,lattice_index(3)) + &
            lattice_vectors(3,lattice_index(2)) * lattice_vectors(3,lattice_index(3))
        scalar_product(3) = &
            lattice_vectors(1,lattice_index(3)) * lattice_vectors(1,lattice_index(1)) + &
            lattice_vectors(2,lattice_index(3)) * lattice_vectors(2,lattice_index(1)) + &
            lattice_vectors(3,lattice_index(3)) * lattice_vectors(3,lattice_index(1))

        if(.not. done) then
            if(maxloc(scalar_product,dim=1) .eq. 1) then
                lattice_vectors(1:3,lattice_index(2)) = &
                    lattice_vectors(1:3,lattice_index(2)) - &
                        lattice_vectors(1:3,lattice_index(1))
            elseif(maxloc(scalar_product,dim=1) .eq. 2) then
                lattice_vectors(1:3,lattice_index(3)) = &
                    lattice_vectors(1:3,lattice_index(3)) - &
                        lattice_vectors(1:3,lattice_index(2))
            elseif(maxloc(scalar_product,dim=1) .eq. 3) then
                lattice_vectors(1:3,lattice_index(3)) = &
                    lattice_vectors(1:3,lattice_index(3)) - &
                        lattice_vectors(1:3,lattice_index(1))
            endif
        endif

        if(done) return
    enddo

    if(.not. done) then
        write(6,*) "Failed delaunay reduction; Stopping program"
        stop
    endif

    ! Adjust lattice vectors and atomic coordinates accordingly
    temp_lattice_vectors = lattice_vectors
    do i = 1, 3
        lattice_vectors(1:3,i) = temp_lattice_vectors(1:3,lattice_index(i))
    enddo
    allocate(temp_atomic_coordinates(3,n_atom))
    temp_atomic_coordinates = atomic_coordinates
    atomic_coordinates(1,1:n_atom) = temp_atomic_coordinates(lattice_index(1),1:n_atom)
    atomic_coordinates(2,1:n_atom) = temp_atomic_coordinates(lattice_index(2),1:n_atom)
    atomic_coordinates(3,1:n_atom) = temp_atomic_coordinates(lattice_index(3),1:n_atom)
    deallocate(temp_atomic_coordinates)

end subroutine delaunayReduction

end module structureinit_m