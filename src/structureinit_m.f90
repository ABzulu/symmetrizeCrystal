module structureinit_m
    use constants, only: eps16, eps6
    use calculatereciprocallattice_m, only: calculateReciprocalLattice

    implicit none

    public :: calculateMetric, restricAtomicCoordinates

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
    if(debug) then
        write(6,*) "restricAtomicCoordinates: lattice_vectors"
        write(6,'(3f16.9)') lattice_vectors(1,1:3)
        write(6,'(3f16.9)') lattice_vectors(2,1:3)
        write(6,'(3f16.9)') lattice_vectors(3,1:3)
        write(6,*) "restricAtomicCoordinates: atomic_coordinates"
        do ia = 1, n_atom
            write(6,'(3f16.9)') atomic_coordinates(1:3,ia)
        enddo
    endif
    call calculateReciprocalLattice(lattice_vectors, reciprocal_lattice_vectors, 0)
    if(debug) then
        write(6,*) "restricAtomicCoordinates: reciprocal_lattice_vectors"
        write(6,'(3f16.9)') reciprocal_lattice_vectors(1,1:3)
        write(6,'(3f16.9)') reciprocal_lattice_vectors(2,1:3)
        write(6,'(3f16.9)') reciprocal_lattice_vectors(3,1:3)
    endif
    if(debug) write(6,*) "restricAtomicCoordinates: atomic_coordinates"
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

subroutine calculateMetric( &
    lattice_vectors, lattice_index, reduced_lattice_vectors, G, debug &
)
    double precision, intent(in) :: lattice_vectors(3,3)
    integer, intent(out) :: lattice_index(3)
    double precision, intent(out) :: reduced_lattice_vectors(3,3), G(3,3)
    logical, intent(in) :: debug

    integer :: temp_i
    double precision :: a_sqr(3), temp_a, temp_lattice_vectors(3,3), scalar_product(3)
    logical :: done

    integer :: ir, i, j

    ! First use delaunay reduction
    temp_lattice_vectors = lattice_vectors
    done = .false.
    lattice_index(1) = 1
    lattice_index(2) = 2
    lattice_index(3) = 3
    do ir = 1, 100
        ! Order lattice vectors so that a_1 <= a_2 <= a_3
        do i = 1, 3
            a_sqr(i) = &
                temp_lattice_vectors(lattice_index(i),1)*temp_lattice_vectors(lattice_index(i),1) + &
                temp_lattice_vectors(lattice_index(i),2)*temp_lattice_vectors(lattice_index(i),2) + &
                temp_lattice_vectors(lattice_index(i),3)*temp_lattice_vectors(lattice_index(i),3)
        enddo
        if(abs(a_sqr(2) - a_sqr(3)) .gt. eps6) then
            temp_a = a_sqr(2)
            a_sqr(2) = a_sqr(3)
            a_sqr(3) = temp_a
            temp_i = lattice_index(2)
            lattice_index(2) = lattice_index(3)
            lattice_index(3) = temp_i
        endif
        if(abs(a_sqr(1) - a_sqr(2)) .gt. eps6) then
            temp_a = a_sqr(1)
            a_sqr(1) = a_sqr(2)
            a_sqr(2) = temp_a
            temp_i = lattice_index(1)
            lattice_index(1) = lattice_index(2)
            lattice_index(2) = temp_i
        endif
        if(abs(a_sqr(2) - a_sqr(3)) .gt. eps6) then
            temp_a = a_sqr(2)
            a_sqr(2) = a_sqr(3)
            a_sqr(3) = temp_a
            temp_i = lattice_index(2)
            lattice_index(2) = lattice_index(3)
            lattice_index(3) = temp_i
        endif

        ! Make sure scalar product of lattice vectors are <= 0
        scalar_product(1) = &
            temp_lattice_vectors(lattice_index(1),1)*temp_lattice_vectors(lattice_index(2),1) + &
            temp_lattice_vectors(lattice_index(1),2)*temp_lattice_vectors(lattice_index(2),2) + &
            temp_lattice_vectors(lattice_index(1),3)*temp_lattice_vectors(lattice_index(2),3)
        scalar_product(2) = &
            temp_lattice_vectors(lattice_index(2),1)*temp_lattice_vectors(lattice_index(3),1) + &
            temp_lattice_vectors(lattice_index(2),2)*temp_lattice_vectors(lattice_index(3),2) + &
            temp_lattice_vectors(lattice_index(2),3)*temp_lattice_vectors(lattice_index(3),3)
        scalar_product(3) = &
            temp_lattice_vectors(lattice_index(3),1)*temp_lattice_vectors(lattice_index(1),1) + &
            temp_lattice_vectors(lattice_index(3),2)*temp_lattice_vectors(lattice_index(1),2) + &
            temp_lattice_vectors(lattice_index(3),3)*temp_lattice_vectors(lattice_index(1),3)

        if(debug) write(6,*) "calculateMetric: scalar product"
        if(debug) write(6,'(3f16.9)') scalar_product(:)

        if((scalar_product(1) .le. eps16) .and. &
           (scalar_product(2) .le. eps16) .and. &
           (scalar_product(3) .le. eps16)) done = .true.

        if(debug) write(6,*) "calculateMetric: temp lattice vectors"
        if(debug) write(6,'(3f16.9)') temp_lattice_vectors(lattice_index(1),1:3)
        if(debug) write(6,'(3f16.9)') temp_lattice_vectors(lattice_index(2),1:3)
        if(debug) write(6,'(3f16.9)') temp_lattice_vectors(lattice_index(3),1:3)

        if(debug) write(6,*) "calculateMetric: done =", done

        if(.not. done) then
            if(maxloc(scalar_product,dim=1) .eq. 1) then
                temp_lattice_vectors(lattice_index(2),1:3) = &
                    temp_lattice_vectors(lattice_index(2),1:3) - &
                    temp_lattice_vectors(lattice_index(1),1:3)
            elseif(maxloc(scalar_product,dim=1) .eq. 2) then
                temp_lattice_vectors(lattice_index(3),1:3) = &
                    temp_lattice_vectors(lattice_index(3),1:3) - &
                    temp_lattice_vectors(lattice_index(2),1:3)
            elseif(maxloc(scalar_product,dim=1) .eq. 3) then
                temp_lattice_vectors(lattice_index(3),1:3) = &
                    temp_lattice_vectors(lattice_index(3),1:3) - &
                    temp_lattice_vectors(lattice_index(1),1:3)
            endif
        endif

        if(done) exit
    enddo

    if(.not. done) then
        write(6,*) "Failed delaunay reduction; Stopping program"
        stop
    endif

    ! Adjust lattice vectors accordingly
    do i = 1, 3
        reduced_lattice_vectors(i,1:3) = temp_lattice_vectors(lattice_index(i),1:3)
    enddo

    ! Calculate metric G = B^TB
    temp_lattice_vectors = G
    do i = 1, 3
        do j = 1, 3
            G(i,j) = &
                reduced_lattice_vectors(i,1) * reduced_lattice_vectors(j,1) + &
                reduced_lattice_vectors(i,2) * reduced_lattice_vectors(j,2) + &
                reduced_lattice_vectors(i,3) * reduced_lattice_vectors(j,3)
        enddo
    enddo
    if(debug) then
        write(6,*) "calculateMetric: G"
        write(6,'(3f16.9)') G(1,1:3)
        write(6,'(3f16.9)') G(2,1:3)
        write(6,'(3f16.9)') G(3,1:3)
    endif

end subroutine calculateMetric

end module structureinit_m