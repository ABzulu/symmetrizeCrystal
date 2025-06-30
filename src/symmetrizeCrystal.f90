program symmetrizeCrystal
    use inlineoptions_m, only: getInlineOptions
    use input_m, only: readInput, formatInput
    use cellReduction_m, only: cellReduction
    use findRotationalSymmetry_m, only: findRotationalSymmetry
    use identifyCrystal_m, only: identifyCrystal
    use restricAtomicCoordinates_m, only: restricAtomicCoordinates
    use findSymmetryOperators_m, only: findSymmetryOperators
    use spacegroup_db, only: loadInSpaceGroup
    use identifySpaceGroup_m, only: identifySpaceGroup
    use symmetrizeStructure_m, only: symmetrizeStructure
    use writeOutput_m, only: writeOutput

    implicit none

    integer :: count_start, count_end, rate
    real :: elapsed_time

    character(132) :: &
        atomic_coordinates_format, input_filename, output_filename
    integer :: n_atom, n_species
    integer, allocatable :: atomic_species_index(:)
    double precision :: &
        lattice_constant, lattice_vectors(3,3), lattice_tol, atomic_tol(3)
    double precision, allocatable :: atomic_coordinates(:,:)
    logical :: debug

    character(132) :: crystal_type
    integer :: W(3,3,48), n_W
    double precision :: reduced_lattice_vectors(3,3)

    double precision :: original_atomic_tol
    integer :: &
        n_symmetry_operators, ix, &
        n_ops(530), rotations(3,3,192,530), translations(3,192,530), ITA_index(530)
    integer, allocatable :: symmetry_operators(:,:,:)
    character(len=17) :: hall_symbol(530)
    logical :: found_space_group

    call system_clock(count_start, rate)

    ! Read in command line options and inputs
    call getInlineOptions( &
        lattice_tol, atomic_tol(1), input_filename, output_filename, debug &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "getInlineOptions took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    ! Read in input.fdf for structure parameters
    call readInput( &
        input_filename, lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates,&
        n_species, atomic_species_index &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "readInput took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call formatInput( &
        lattice_constant, lattice_vectors, &
        atomic_coordinates_format, atomic_coordinates, n_atom, &
        debug &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "formatInput took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call cellReduction( &
        lattice_vectors, reduced_lattice_vectors, lattice_tol, debug &
    )
    
    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "cellReduction took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call findRotationalSymmetry(reduced_lattice_vectors, lattice_tol, W, n_W, debug)

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "findRotationalSymmetry took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call identifyCrystal(W, n_W, crystal_type, debug)

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "identifyCrystal took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call restricAtomicCoordinates( &
        reduced_lattice_vectors, n_atom, atomic_coordinates, debug &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "restricAtomicCoordinates took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    original_atomic_tol = atomic_tol(1)
    do ix = 1, 3
        atomic_tol(ix) = &
            original_atomic_tol / &
            sqrt( &
                reduced_lattice_vectors(ix,1)*reduced_lattice_vectors(ix,1) + &
                reduced_lattice_vectors(ix,2)*reduced_lattice_vectors(ix,2) + &
                reduced_lattice_vectors(ix,3)*reduced_lattice_vectors(ix,3) &
            )
        if(debug) &
            write(6,'(a,i4,a,f20.14)') "Atomic tolerance in ", ix, " axis =", atomic_tol(ix)
    enddo

    call findSymmetryOperators( &
        reduced_lattice_vectors, &
        n_atom, atomic_coordinates, n_species, atomic_species_index, & 
        atomic_tol, n_W, W, n_symmetry_operators, symmetry_operators, debug &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "findSymmetryOperators took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call loadInSpaceGroup( &
        n_ops, rotations, translations, hall_symbol, ITA_index &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "loadInSpaceGroup took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call identifySpaceGroup( &
        symmetry_operators, n_symmetry_operators, &
        n_ops, rotations, translations, hall_symbol, ITA_index, found_space_group, debug &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "identifySpaceGroup took ", elapsed_time, " seconds"

    if(found_space_group) then
        continue
    else
        write(6,'(a)') "Could not find space group."
        stop
    endif
            
    call system_clock(count_start, rate)

    call symmetrizeStructure( &
        reduced_lattice_vectors, n_atom, atomic_coordinates, atomic_species_index, &
        n_symmetry_operators, symmetry_operators, atomic_tol, debug &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "symmetrizeStructure took ", elapsed_time, " seconds"

    call system_clock(count_start, rate)

    call writeOutput( &
        lattice_constant, lattice_vectors, &
        n_atom, atomic_coordinates_format, atomic_coordinates, atomic_species_index, &
        reduced_lattice_vectors, output_filename &
    )

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') "writeOutput took ", elapsed_time, " seconds"

    ! Deallocate all allocatable arrays
    deallocate(atomic_species_index, atomic_coordinates)

end program symmetrizeCrystal