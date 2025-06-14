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

    integer :: iounit, i, ia, ix, jx, T(3,3)
    double precision :: temp_a(3), reciprocal_lattice_vector(3,3), temp_lattice_vector(3,3), temp_TT(3,3)
    double precision, allocatable :: temp_atomic_coordinates(:,:)
    logical :: leqi

    open(newunit=iounit, file=output_filename, status='replace')

    write(iounit,'(a,f20.10,a)') "LatticeConstant", lattice_constant, " Bohr"

    ! Change atomic coordinates to the original format
    allocate(temp_atomic_coordinates(3,n_atom))
    temp_atomic_coordinates = atomic_coordinates

    write(6,'(a)') "writeOutput: Reduced lattice vectors"
    write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
    write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
    write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
    write(6,'(a)') "writeOutput: Atomic coordinates in fractional coordinates"
    do ia = 1, n_atom
        write(6,'(3f16.9)') temp_atomic_coordinates(1:3,ia)
    enddo

    call calculateReciprocalLattice(lattice_vector, reciprocal_lattice_vector, 0)

    do ix = 1, 3;do jx = 1, 3
        T(ix,jx) = nint( &
            reciprocal_lattice_vector(ix,1) * reduced_lattice_vectors(jx,1) + &
            reciprocal_lattice_vector(ix,2) * reduced_lattice_vectors(jx,2) + &
            reciprocal_lattice_vector(ix,3) * reduced_lattice_vectors(jx,3) &
        )
    enddo;enddo
    write(6,'(a)') "writeOutput: T"
    write(6,'(3i4)') T(1,1:3)
    write(6,'(3i4)') T(2,1:3)
    write(6,'(3i4)') T(3,1:3)
    do ia = 1, n_atom
        temp_a(1:3) = temp_atomic_coordinates(1:3,ia)
        do ix = 1, 3
            temp_atomic_coordinates(ix,ia) = &
                T(ix,1) * temp_a(1) + &
                T(ix,2) * temp_a(2) + &
                T(ix,3) * temp_a(3)
            temp_atomic_coordinates(ix,ia) = &
                modulo(temp_atomic_coordinates(ix,ia),1.d0)
        enddo
    enddo

    call calculateReciprocalLattice(reduced_lattice_vectors, reciprocal_lattice_vector, 0)

    do ix = 1, 3;do jx = 1, 3
        T(ix,jx) = nint( &
            reciprocal_lattice_vector(ix,1) * lattice_vector(jx,1) + &
            reciprocal_lattice_vector(ix,2) * lattice_vector(jx,2) + &
            reciprocal_lattice_vector(ix,3) * lattice_vector(jx,3) &
        )
    enddo;enddo
    write(6,'(a)') "writeOutput: lattice_vector"
    write(6,'(3f16.9)') lattice_vector(1,1:3)
    write(6,'(3f16.9)') lattice_vector(2,1:3)
    write(6,'(3f16.9)') lattice_vector(3,1:3)
    write(6,'(a)') "writeOutput: reduced_lattice_vectors"
    write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
    write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
    write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
    write(6,'(a)') "writeOutput: T"
    write(6,'(3i4)') T(1,1:3)
    write(6,'(3i4)') T(2,1:3)
    write(6,'(3i4)') T(3,1:3)
    do ix = 1, 3
        do jx = 1, 3
            temp_lattice_vector(ix,jx) = &
                reduced_lattice_vectors(1,ix)*T(1,jx) + &
                reduced_lattice_vectors(2,ix)*T(2,jx) + &
                reduced_lattice_vectors(3,ix)*T(3,jx)
        enddo
    enddo

    write(iounit,'(a)') "%block    LatticeVectors"
    do i = 1,3
        write(iounit,'(3f16.10)') temp_lattice_vector(i,1:3)/lattice_constant
    enddo
    write(iounit,'(a)') "%endlock  LatticeVectors"

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
                    lattice_vector(1,i) * temp_a(1) + &
                    lattice_vector(2,i) * temp_a(2) + &
                    lattice_vector(3,i) * temp_a(3)
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

    write(iounit,'(2a)') "AtomicCoordinatesFormat        ", atomic_coordinates_format
    write(iounit,'(a)') "%block    AtomicCoordinatesAndAtomicSpecies"
    do ia = 1, n_atom
        write(iounit,'(3f16.10,i4)') &
            temp_atomic_coordinates(1:3,ia), atomic_species_index(ia)
    enddo
    write(iounit,'(a)') "%endblock AtomicCoordinatesAndAtomicSpecies"

    deallocate(temp_atomic_coordinates)

end subroutine

end module writeOutput_m