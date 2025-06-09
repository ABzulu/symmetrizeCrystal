program symmetrizeCrystal
    use inlineoptions_m, only: getInlineOptions
    use input_m, only: readInput, formatInput
    use structureinit_m, only: restricAtomicCoordinates, calculateMetric
    use symmetry_m, only: findRotationalSymmetry, findTranslationalSymmetry
    use output_m, only: writeOutput

    implicit none

    logical :: debug

    character(132) :: atomic_coordinates_format, output_filename
    integer :: n_atom
    integer, allocatable :: atomic_species_index(:)
    double precision :: &
        lattice_constant, lattice_vectors(3,3), lattice_tol, atomic_tol
    double precision, allocatable :: atomic_coordinates(:,:)

    integer :: lattice_index(3), n_W, n_symm_op
    double precision :: reduced_lattice_vectors(3,3), G(3,3)
    double precision, allocatable :: symm_op(:,:,:)
    integer :: W(3,3,192) ! Rotational matrices

    ! Read in command line options and inputs
    call getInlineOptions(lattice_tol, atomic_tol, output_filename, debug)

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

    if(debug) write(6,*) "Checkpoint: Read input"

    call calculateMetric( &
        lattice_vectors, lattice_index, reduced_lattice_vectors, G, debug &
    )
    call restricAtomicCoordinates(reduced_lattice_vectors, n_atom, atomic_coordinates, debug)
    call findRotationalSymmetry(G, lattice_tol, W, n_W, debug)
    ! Atomic coordinates are changed to fractional between -0.5 < x <= 0.5
    call findTranslationalSymmetry( &
        n_atom, atomic_coordinates, atomic_species_index, atomic_tol,&
        n_W, W, n_symm_op, symm_op, debug &
    )

    if(debug) write(6,*) "Checkpoint: Find symmetry"

    call writeOutput( &
        lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
        output_filename &
    )

    ! Deallocate all allocatable arrays
    deallocate(atomic_species_index, atomic_coordinates)

end program symmetrizeCrystal