module output_m
    use constants, only: ANGSTROM_AU
    use calculatereciprocallattice_m, only: calculateReciprocalLattice

    implicit none

    public :: writeOutput

    private

contains

subroutine writeOutput( &
    lattice_constant, lattice_vector, &
    n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
    output_filename &
)
    integer, intent(in) :: n_atom, atomic_species_index(n_atom)
    double precision, intent(in) :: lattice_constant, lattice_vector(3,3)
    double precision, intent(in) :: atomic_coordinates(3,n_atom)
    character(len=132), intent(in) :: atomic_coordinates_format, output_filename

    integer :: iounit, i, ia
    double precision :: reciprocal_lattice_vector(3,3)
    double precision, allocatable :: temp_atomic_coordinates(:,:)
    logical :: leqi

    open(newunit=iounit, file=output_filename, status='replace')

    write(iounit,'(a,f20.10)') "LatticeConstant", lattice_constant

    write(iounit,'(a)') "%block    LatticeVectors"
    do i = 1,3
        write(iounit,'(3f16.10)') lattice_vector(1:3,i)
    enddo
    write(iounit,'(a)') "%endlock  LatticeVectors"

    ! Change atomic coordinates to the original format
    allocate(temp_atomic_coordinates(3,n_atom))
    temp_atomic_coordinates = atomic_coordinates

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
    elseif(leqi(atomic_coordinates_format,'ScaledByLatticeVectors') .or. &
           leqi(atomic_coordinates_format,'Fractional')) then
        call calculateReciprocalLattice(lattice_vector, reciprocal_lattice_vector, 0)
        do ia = 1, n_atom
            do i = 1, 3
                temp_atomic_coordinates(i,ia) = &
                    reciprocal_lattice_vector(1,i) * atomic_coordinates(1,ia) + &
                    reciprocal_lattice_vector(2,i) * atomic_coordinates(2,ia) + &
                    reciprocal_lattice_vector(3,i) * atomic_coordinates(3,ia)
            enddo
        enddo
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

end module output_m