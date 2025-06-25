module findSymmetryOperators_m
    use constants, only: eps16
    use findDuplicateSymmetryOperators_m, only: findDuplicateSymmetryOperators

    implicit none

    public :: findSymmetryOperators

    private

contains

subroutine findSymmetryOperators( &
    lattice_vectors, n_atom, atomic_coordinates, n_species, atomic_species_index, &
    atomic_tol, n_W, W, n_symmetry_operators, symmetry_operators, debug &
)
    integer, intent(in)               :: n_atom, n_species, atomic_species_index(n_atom)
    double precision, intent(in)      :: &
        lattice_vectors(3,3), atomic_coordinates(3,n_atom), atomic_tol(3)
    logical, intent(in)               :: debug
    integer, intent(inout)            :: n_W, W(3,3,48)
    integer, intent(out)              :: n_symmetry_operators
    integer, allocatable, intent(out) :: symmetry_operators(:,:,:)

    integer              :: is, ia, n_candidates, counter, iw, ja, ix, ic, jx, kx, iix
    integer, allocatable :: &
        n_atom_per_species(:), candidate_symmetry_operators(:,:,:), atom_index(:)
    double precision     :: t
    logical, allocatable :: is_symmetric(:)

    integer                       :: iy, iz
    double precision, allocatable :: images(:,:)
    logical                       :: found_match

    integer :: count_start, count_end, rate
    real    :: elapsed_time
    
    call system_clock(count_start, rate)

    allocate(n_atom_per_species(n_species))
    n_atom_per_species(:) = 0
    do is = 1, n_species
        do ia = 1, n_atom
            if(atomic_species_index(ia) .eq. is) &
                n_atom_per_species(is) = n_atom_per_species(is) + 1
        enddo
    enddo

    n_candidates = 0
    do is = 1, n_species
        n_candidates = &
            n_candidates + n_atom_per_species(is)*n_atom_per_species(is)
    enddo
    if(n_candidates .lt. 1728) then
        n_candidates = n_candidates * n_W
        allocate(candidate_symmetry_operators(3,4,n_candidates))

        ! Find candidate translational operator from other atoms of the same species
        counter = 0
        if(debug) write(6,'(a)') "findSymmetryOperators: candidate_translational_operators"
        do iw = 1, n_W
            do ia = 1, n_atom
                do ja = 1, n_atom
                    if(.not. (atomic_species_index(ia) .eq. atomic_species_index(ja))) cycle
                    counter = counter + 1
                    candidate_symmetry_operators(1:3,1:3,counter) = W(1:3,1:3,iw)

                    do ix = 1, 3
                        t = &
                            atomic_coordinates(ix,ja) - ( &
                                dble(W(ix,1,iw)) * atomic_coordinates(1,ia) + &
                                dble(W(ix,2,iw)) * atomic_coordinates(2,ia) + &
                                dble(W(ix,3,iw)) * atomic_coordinates(3,ia) &                            
                            )
                        t = modulo(t,1.d0)
                        candidate_symmetry_operators(ix,4,counter) = nint(t * 12.d0)
                        if(candidate_symmetry_operators(ix,4,counter) .eq. 12) &
                            candidate_symmetry_operators(ix,4,counter) = 0
                    enddo
                    if(debug) then
                        write(6,'(a,i8)') "findSymmetryOperators: ic = ", counter
                        write(6,'(4i4)') candidate_symmetry_operators(1,1:4,counter)
                        write(6,'(4i4)') candidate_symmetry_operators(2,1:4,counter)
                        write(6,'(4i4)') candidate_symmetry_operators(3,1:4,counter)
                    endif
                enddo
            enddo
        enddo
        
        call findDuplicateSymmetryOperators( &
            lattice_vectors, candidate_symmetry_operators, n_candidates, debug &
        )

        allocate(is_symmetric(n_candidates))

        if(debug) then
            write(6,'(a,i4)') "findSymmetryOperators: Number of candidates =", n_candidates
            do ic = 1, n_candidates
                write(6,'(a,i4)') "findSymmetryOperators: Candidate index =", ic
                write(6,'(4i4)') candidate_symmetry_operators(1,1:4,ic)
                write(6,'(4i4)') candidate_symmetry_operators(2,1:4,ic)
                write(6,'(4i4)') candidate_symmetry_operators(3,1:4,ic)
            enddo
        endif
    else
        n_candidates = 1728 * n_W
        allocate(candidate_symmetry_operators(3,4,n_candidates))
        counter = 0
        do iw = 1, n_W
            do ix = 0, 11;do jx = 0, 11;do kx = 0, 11
                counter = counter + 1
                do iix = 1, 3
                    candidate_symmetry_operators(iix,1:3,counter) = W(iix,1:3,iw)
                enddo
                candidate_symmetry_operators(1,4,counter) = ix
                candidate_symmetry_operators(2,4,counter) = jx
                candidate_symmetry_operators(3,4,counter) = kx
            enddo;enddo;enddo
        enddo

        ! call findDuplicateSymmetryOperators( &
        !     lattice_vectors, candidate_symmetry_operators, n_candidates, debug &
        ! )        

        allocate(is_symmetric(n_candidates))

        if(debug) write(6,'(a)') "findSymmetryOperators: Searching through all possible Symmetry operators"
    endif

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') &
        "findSymmetryOperators-make candidate symmetry operators took ", &
        elapsed_time, " seconds"

    call system_clock(count_start, rate)

    allocate(images(3,n_atom*125), atom_index(n_atom*125))
    ! Make image using symmetry operator
    n_symmetry_operators = 0
    is_symmetric(1:n_candidates) = .false.
    do ic = 1, n_candidates
        do ia = 1, n_atom
            atom_index(ia) = ia
            do ix = 1, 3
                images(ix,ia) = &
                    dble(candidate_symmetry_operators(ix,1,ic)) * atomic_coordinates(1,ia) + &
                    dble(candidate_symmetry_operators(ix,2,ic)) * atomic_coordinates(2,ia) + &
                    dble(candidate_symmetry_operators(ix,3,ic)) * atomic_coordinates(3,ia) + &
                    dble(candidate_symmetry_operators(ix,4,ic))/12.d0
            enddo
        enddo

        ! Make supercell from the image
        counter = n_atom
        do ix = -2, 2;do iy = -2, 2;do iz = -2, 2
            if((ix .eq. 0) .and. (iy .eq. 0) .and. (iz .eq. 0)) cycle

            do ia = 1, n_atom
                counter = counter + 1
                atom_index(counter) = ia
                images(1,counter) = images(1,ia) + ix
                images(2,counter) = images(2,ia) + iy
                images(3,counter) = images(3,ia) + iz
            enddo
        enddo;enddo;enddo

        ! Check if all atoms have matching image
        do ia = 1, n_atom
            found_match = .false.
            do ja = 1, n_atom*125
                if(atomic_species_index(ia) .ne. atomic_species_index(atom_index(ja))) cycle
                if( &
                    (abs(atomic_coordinates(1,ia) - images(1,ja)) .lt. atomic_tol(1)) .and. &
                    (abs(atomic_coordinates(2,ia) - images(2,ja)) .lt. atomic_tol(2)) .and. &
                    (abs(atomic_coordinates(3,ia) - images(3,ja)) .lt. atomic_tol(3)) &
                ) then
                    found_match = .true.
                    exit
                endif
            enddo
            if(.not. found_match) exit
        enddo

        if(.not. found_match) cycle

        n_symmetry_operators = n_symmetry_operators + 1
        is_symmetric(ic) = .true.
    enddo

    allocate(symmetry_operators(3,4,n_symmetry_operators))

    write(6,'(a,i6)') "findSymmetryOperators: Number of Symmetry operators =", n_symmetry_operators
    write(6,'(a)') "findSymmetryOperators: symmetry_operators"
    counter = 0
    do ic = 1, n_candidates
        if(is_symmetric(ic)) then
            counter = counter + 1
            symmetry_operators(1:3,1:4,counter) = candidate_symmetry_operators(1:3,1:4,ic)
            write(6,'(a,i6)') "symmetry operator index", counter
            write(6,'(4i4)') symmetry_operators(1,1:4,counter)
            write(6,'(4i4)') symmetry_operators(2,1:4,counter)
            write(6,'(4i4)') symmetry_operators(3,1:4,counter)
        endif
    enddo

    call system_clock(count_end)
    elapsed_time = real(count_end - count_start) / real(rate)
    write(6,'(a,f14.4,a)') &
        "findSymmetryOperators-find symmetry operators using images took ", &
        elapsed_time, " seconds"

    deallocate(n_atom_per_species)
    deallocate(candidate_symmetry_operators, is_symmetric, images)

end subroutine findSymmetryOperators

end module findSymmetryOperators_m