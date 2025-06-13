module writeOutput_m
    use constants, only: ANGSTROM_AU, eps16
    use calculatereciprocallattice_m, only: calculateReciprocalLattice

    implicit none

    public :: writeOutput

    private

contains

subroutine writeOutput( &
    lattice_constant, lattice_vector, &
    n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
    reduced_lattice_vectors, output_filename &
)
    integer, intent(in) :: n_atom, atomic_species_index(n_atom)
    double precision, intent(in) :: &
        lattice_constant, lattice_vector(3,3), reduced_lattice_vectors(3,3), &
        atomic_coordinates(3,n_atom)
    character(len=132), intent(in) :: atomic_coordinates_format, output_filename

    integer :: iounit, i, ia
    double precision :: temp_a(3), reciprocal_lattice_vector(3,3)
    double precision, allocatable :: temp_atomic_coordinates(:,:)
    logical :: leqi

    open(newunit=iounit, file=output_filename, status='replace')

    write(iounit,'(a,f20.10,a)') "LatticeConstant", lattice_constant, " Bohr"

    write(iounit,'(a)') "%block    LatticeVectors"
    do i = 1,3
        write(iounit,'(3f16.10)') lattice_vector(i,1:3)/lattice_constant
    enddo
    write(iounit,'(a)') "%endlock  LatticeVectors"

    ! Change atomic coordinates to the original format
    allocate(temp_atomic_coordinates(3,n_atom))
    temp_atomic_coordinates = atomic_coordinates

    temp_atomic_coordinates(:,:) = temp_atomic_coordinates(:,:) + 0.5d0

    write(6,'(a)') "writeOutput: Reduced lattice vectors"
    write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
    write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
    write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
    write(6,'(a)') "writeOutput: Atomic coordinates in fractional coordinates"
    do ia = 1, n_atom
        write(6,'(3f16.9)') temp_atomic_coordinates(1:3,ia)
    enddo

    ! Change the coordinates in Bohr
    do ia = 1, n_atom
        temp_a(1:3) = temp_atomic_coordinates(1:3,ia)
        do i = 1, 3
            temp_atomic_coordinates(i,ia) = &
                reduced_lattice_vectors(i,1) * temp_a(1) + &
                reduced_lattice_vectors(i,2) * temp_a(2) + &
                reduced_lattice_vectors(i,3) * temp_a(3)
        enddo
    enddo

    ! Change the coordinates in Fractional w.r.t input lattice vectors
    call calculateReciprocalLattice( &
        lattice_vector, reciprocal_lattice_vector, 0 &
    )
    do ia = 1, n_atom
        temp_a(1:3) = temp_atomic_coordinates(1:3,ia)
        do i = 1, 3
            temp_atomic_coordinates(i,ia) = &
                reciprocal_lattice_vector(i,1) * temp_a(1) + &
                reciprocal_lattice_vector(i,2) * temp_a(2) + &
                reciprocal_lattice_vector(i,3) * temp_a(3)
            temp_atomic_coordinates(i,ia) = &
                modulo(temp_atomic_coordinates(i,ia),1.d0)
            ! This is to take care of exactly 1.0d0 
            if( temp_atomic_coordinates(i,ia) - 1.0d0 .gt. eps16) then
                temp_atomic_coordinates(i,ia) = 0.0d0
            endif
        enddo
    enddo

    write(6,'(a)') "writeOutput: Lattice vectors"
    write(6,'(3f16.9)') lattice_vector(1,1:3)
    write(6,'(3f16.9)') lattice_vector(2,1:3)
    write(6,'(3f16.9)') lattice_vector(3,1:3)
    write(6,'(a)') "writeOutput: Atomic coordinates in fractional coordinates"
    do ia = 1, n_atom
        write(6,'(3f16.9)') temp_atomic_coordinates(1:3,ia)
    enddo

    if(leqi(atomic_coordinates_format,'ScaledByLatticeVectors') .or. &
       leqi(atomic_coordinates_format,'Fractional')) then
        continue
    else
        do ia = 1, n_atom
            temp_a(1:3) = temp_atomic_coordinates(1:3,ia)
            do i = 1, 3
                temp_atomic_coordinates(i,ia) = &
                    lattice_vector(i,1) * temp_a(1) + &
                    lattice_vector(i,2) * temp_a(2) + &
                    lattice_vector(i,3) * temp_a(3)
            enddo
        enddo
        if(leqi(atomic_coordinates_format,'NotScaledCartesianBohr') .or. &
           leqi(atomic_coordinates_format,'Bohr')) then
            continue
        elseif(leqi(atomic_coordinates_format,'NotScaledCartesianAng') .or. &
               leqi(atomic_coordinates_format,'Ang')) then
            temp_atomic_coordinates = temp_atomic_coordinates / ANGSTROM_AU
        elseif(leqi(atomic_coordinates_format,'ScaledCartesian')) then
            if(lattice_constant .eq. 0.d0) then
                write(6,*) "Lattice constant not given when atomic coordinates format is 'ScaledCartesian'; Stoping program"
                stop
            endif
            temp_atomic_coordinates = temp_atomic_coordinates / lattice_constant
        endif
    endif

    write(iounit,'(2a)') "AtomicCoordinatesFormat", atomic_coordinates_format
    write(iounit,'(a)') "%block    AtomicCoordinatesAndAtomicSpecies"
    do ia = 1, n_atom
        write(iounit,'(3f16.10,i4)') &
            temp_atomic_coordinates(1:3,ia), atomic_species_index(ia)
    enddo
    write(iounit,'(a)') "%endblock AtomicCoordinatesAndAtomicSpecies"

    deallocate(temp_atomic_coordinates)

end subroutine

end module writeOutput_m