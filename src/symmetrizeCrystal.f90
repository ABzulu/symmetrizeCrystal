program symmetrizeCrystal
    use inlineoptions_m, only: getInlineOptions
    use input_m, only: readInput, formatInput
    use output_m, only: writeOutput

    implicit none

    logical :: debug

    character(132) :: atomic_coordinates_format, output_filename
    integer :: n_atom
    integer, allocatable :: atomic_species_index(:)
    double precision :: lattice_constant, lattice_vector(3,3)
    double precision, allocatable :: atomic_coordinates(:,:)

    ! Read in command line options and inputs
    call getInlineOptions(output_filename, debug)

    ! Read in input.fdf for structure parameters
    call readInput( &
        lattice_constant, lattice_vector, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index &
    )
    ! Change atomic coordinates to Bohr
    call formatInput( &
        lattice_constant, lattice_vector, &
        atomic_coordinates_format, atomic_coordinates, n_atom, debug &
    )

    call writeOutput( &
        lattice_constant, lattice_vector, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
        output_filename &
    )

    ! Deallocate all allocatable arrays
    deallocate(atomic_species_index, atomic_coordinates)

end program symmetrizeCrystal