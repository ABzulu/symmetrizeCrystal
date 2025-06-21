program symmetrizeCrystal
    use inlineoptions_m, only: getInlineOptions
    use input_m, only: readInput, formatInput
    use niggliReduction_m, only: niggliReduction
    use findRotationalSymmetry_m, only: findRotationalSymmetry
    use identifyCrystal_m, only: identifyCrystal
    use toITAConventionalCell_m, only: toITAConventionalCell
    use restricAtomicCoordinates_m, only: restricAtomicCoordinates
    use findTranslationalSymmetry_m, only: findTranslationalSymmetry
    use spacegroup_db, only: loadInSpaceGroup
    use identifySpaceGroup_m, only: identifySpaceGroup
    use findDuplicateTranslations_m, only: findDuplicateTranslations
    use symmetrizeConventionalCell_m, only: symmetrizeConventionalCell
    use writeOutput_m, only: writeOutput

    implicit none

    logical :: debug

    character(132) :: &
        atomic_coordinates_format, input_filename, output_filename
    integer :: n_atom, n_species, maxIteration
    integer, allocatable :: atomic_species_index(:)
    double precision :: &
        lattice_constant, lattice_vectors(3,3), lattice_tol, atomic_tol(3)
    double precision, allocatable :: atomic_coordinates(:,:)

    character(132) :: crystal_type
    integer :: W(3,3,48), n_W, W_type(48)
    double precision :: reduced_lattice_vectors(3,3)

    double precision :: ITA_lattice_vectors(3,3), original_atomic_tol
    integer :: &
        n_symm_op, ix, &
        n_ops(530), rotations(3,3,192,530), translations(3,192,530), ITA_index(530)
    integer, allocatable :: symm_op(:,:,:)
    character(len=17) :: hall_symbol(530)
    logical :: found_space_group

    integer :: ir

    ! Read in command line options and inputs
    call getInlineOptions( &
        lattice_tol, atomic_tol(1), input_filename, output_filename, &
        maxIteration, debug &
    )

    ! Read in input.fdf for structure parameters
    call readInput( &
        input_filename, lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates,&
        n_species, atomic_species_index &
    )

    call formatInput( &
        lattice_constant, lattice_vectors, &
        atomic_coordinates_format, atomic_coordinates, n_atom, &
        debug &
    )

    if(debug) write(6,'(a)') "Checkpoint: Read input"

    call loadInSpaceGroup( &
        n_ops, rotations, translations, hall_symbol, ITA_index &
    )

    call niggliReduction( &
        lattice_vectors, reduced_lattice_vectors, lattice_tol, debug &
    )

    call findRotationalSymmetry(reduced_lattice_vectors, lattice_tol, W, n_W, debug)

    call identifyCrystal(W, n_W, W_type, crystal_type, debug)

    if(debug) write(6,'(a)') "Checkpoint: Determined point group"

    call toITAConventionalCell( &
        crystal_type, reduced_lattice_vectors, ITA_lattice_vectors, debug &
    )

    call findRotationalSymmetry(ITA_lattice_vectors, lattice_tol, W, n_W, debug)

    call restricAtomicCoordinates( &
        ITA_lattice_vectors, lattice_vectors, n_atom, atomic_coordinates, debug &
    )

    original_atomic_tol = atomic_tol(1)
    do ix = 1, 3
        atomic_tol(ix) = &
            original_atomic_tol / &
            sqrt( &
                ITA_lattice_vectors(ix,1)*ITA_lattice_vectors(ix,1) + &
                ITA_lattice_vectors(ix,2)*ITA_lattice_vectors(ix,2) + &
                ITA_lattice_vectors(ix,3)*ITA_lattice_vectors(ix,3) &
            )
    enddo

    if(debug) write(6,'(a)') "Checkpoint: Made compatible structure"

    call findTranslationalSymmetry( &
        n_atom, atomic_coordinates, n_species, atomic_species_index, & 
        atomic_tol, n_W, W, n_symm_op, symm_op, debug &
    )

    call findDuplicateTranslations( &
        ITA_lattice_vectors, symm_op, n_symm_op, debug &
    )

    call identifySpaceGroup( &
        symm_op, n_symm_op, &
        n_ops, rotations, translations, hall_symbol, found_space_group, debug &
    )

    if(found_space_group) then
        continue
    else
        write(6,'(a)') "Could not find space group."
        stop
    endif
            
    if(debug) write(6,'(a)') "Checkpoint: Determined space group"

    call symmetrizeConventionalCell( &
        ITA_lattice_vectors, n_atom, atomic_coordinates, atomic_species_index, &
        n_symm_op, symm_op, atomic_tol, debug &
    )

    if(debug) write(6,'(a)') "Checkpoint: Symmetrized structure"

    call writeOutput( &
        lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
        ITA_lattice_vectors, output_filename &
    )

    ! Deallocate all allocatable arrays
    deallocate(atomic_species_index, atomic_coordinates)

end program symmetrizeCrystal