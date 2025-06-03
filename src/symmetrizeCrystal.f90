program symmetrizeCrystal
    use inlineoptions_m, only: getInlineOptions
    use input_m, only: readInput, formatInput
    use structureinit_m, only: restricAtomicCoordinates, delaunayReduction
    use output_m, only: writeOutput

    implicit none

    logical :: debug

    character(132) :: atomic_coordinates_format, output_filename
    integer :: n_atom
    integer, allocatable :: atomic_species_index(:)
    double precision :: lattice_constant, lattice_vectors(3,3)
    double precision, allocatable :: atomic_coordinates(:,:)

    integer :: lattice_index(3)

    ! Read in command line options and inputs
    call getInlineOptions(output_filename, debug)

    ! Read in input.fdf for structure parameters
    call readInput( &
        lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index &
    )
    ! Change atomic coordinates to Bohr
    ! This is here because the lattice vectors are automatically transformed to Bohr
    call formatInput( &
        lattice_constant, lattice_vectors, &
        atomic_coordinates_format, atomic_coordinates, n_atom, debug &
    )

    ! Atomic coordinates are changed to fractional between -0.5 < x <= 0.5
    call restricAtomicCoordinates(lattice_vectors, n_atom, atomic_coordinates)
    call delaunayReduction( &
        lattice_vectors, n_atom, atomic_coordinates, lattice_index &
    )

    call writeOutput( &
        lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
        output_filename &
    )

    ! Deallocate all allocatable arrays
    deallocate(atomic_species_index, atomic_coordinates)

end program symmetrizeCrystal