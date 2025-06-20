module findTranslationalSymmetry_m
    use constants, only: eps16

    implicit none

    public :: findTranslationalSymmetry

    private

contains

subroutine findTranslationalSymmetry( &
    n_atom, atomic_coordinates, n_species, atomic_species_index, &
    atomic_tol, n_W, W, n_symm_op, symm_op, debug &
)
    integer, intent(in) :: n_atom, n_species, atomic_species_index(n_atom)
    double precision, intent(in) :: atomic_coordinates(3,n_atom), atomic_tol(3)
    logical, intent(in) :: debug
    integer, intent(inout) :: n_W, W(3,3,48)
    integer, intent(out) :: n_symm_op
    integer, allocatable, intent(out) :: symm_op(:,:,:)

    integer :: is, n_candidates, ia, ja, iw, ix, jx, kx, counter, ic
    integer, allocatable :: n_atom_per_species(:), W_index(:)
    double precision :: delta_x(3)
    double precision, allocatable :: candidate_translational_operator(:,:)
    logical :: found_match, W_match(48)
    logical, allocatable :: is_symmetric(:)

    
    ! When there are more than 500 atoms just use brute force search.
    if(n_atom .lt. 500) then
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
        n_candidates = n_candidates * n_W
        allocate(candidate_translational_operator(3,n_candidates), is_symmetric(n_candidates))
        allocate(W_index(n_candidates))

        ! Find candidate translational operator from other atoms of the same species
        counter = 0
        ! if(debug) write(6,'(a)') "findTranslationalSymmetry: candidate_translational_operators"
        do iw = 1, n_W
            do ia = 1, n_atom
                do ja = 1, n_atom
                    if(.not. (atomic_species_index(ia) .eq. atomic_species_index(ja))) cycle
                    counter = counter + 1
                    W_index(counter) = iw
                    do ix = 1, 3
                        candidate_translational_operator(ix,counter) = &
                            atomic_coordinates(ix,ja) - ( &
                            dble(W(ix,1,iw)) * atomic_coordinates(1,ia) + &
                            dble(W(ix,2,iw)) * atomic_coordinates(2,ia) + &
                            dble(W(ix,3,iw)) * atomic_coordinates(3,ia))
                        candidate_translational_operator(ix,counter) = &
                            modulo(candidate_translational_operator(ix,counter),1.d0)
                        if(abs(candidate_translational_operator(ix,counter) - 1.d0) .lt. eps16) &
                            candidate_translational_operator(ix,counter) = 0.d0
                    enddo
                    ! if(debug) write(6,'(i6,3f16.9)') &
                    !     counter, candidate_translational_operator(1:3,counter)
                enddo
            enddo
        enddo        
    else
        n_candidates = 1728 * n_W
        allocate(candidate_translational_operator(3,n_candidates))
        counter = 0
        do ix = 0, 11;do jx = 0,11;do kx = 0,11
            counter = counter + 1
            candidate_translational_operator(1,counter) = dble(ix) / 12.d0
            candidate_translational_operator(2,counter) = dble(jx) / 12.d0
            candidate_translational_operator(3,counter) = dble(kx) / 12.d0
        enddo;enddo;enddo
    endif

    W_match(:) = .false.
    ! Check if the candidate translational operator maps all the atoms to the
    ! same species of atoms.
    ! if(debug) write(6,*) "findTranslationalSymmetry: delta_x"
    is_symmetric(:) = .false.
    do ic = 1, n_candidates
        if(debug) then
            if( &
                (W(1,1,W_index(ic)) .eq. 1) .and. &
                (W(2,2,W_index(ic)) .eq. 1) .and. &
                (W(3,3,W_index(ic)) .eq. 1) .and. &
                (W(1,2,W_index(ic)) .eq. 0) .and. &
                (W(2,1,W_index(ic)) .eq. 0) .and. &
                (W(1,3,W_index(ic)) .eq. 0) .and. &
                (W(3,1,W_index(ic)) .eq. 0) .and. &
                (W(2,3,W_index(ic)) .eq. 0) .and. &
                (W(3,2,W_index(ic)) .eq. 0) &
            ) then
                ! write(6,'(a,3f16.9)') "Candidate translation operator", candidate_translational_operator(1:3,ic)
            endif
        endif
        do ia = 1, n_atom
            found_match = .false.
            do ja = 1, n_atom
                if(.not. (atomic_species_index(ia) .eq. atomic_species_index(ja))) cycle
                do ix = 1, 3
                    delta_x(ix) = &
                        dble(W(ix,1,W_index(ic))) * atomic_coordinates(1,ia) + &
                        dble(W(ix,2,W_index(ic))) * atomic_coordinates(2,ia) + &
                        dble(W(ix,3,W_index(ic))) * atomic_coordinates(3,ia) + &
                        candidate_translational_operator(ix,ic)
                    delta_x(ix) = modulo(delta_x(ix) + 0.5d0, 1.d0) - 0.5d0
                    if(abs(delta_x(ix) + 0.5d0) .lt. eps16 ) delta_x(ix) = 0.5d0

                    delta_x(ix) = abs(delta_x(ix) - atomic_coordinates(ix,ja))
                enddo

                ! if(debug) write(6,'(a,3i6,3f16.9)') &
                !     "ic, ia, ja, delta_x = ", ic, ia, ja, delta_x(1:3)
                if( &
                    (delta_x(1) .lt. atomic_tol(1)) .and. &
                    (delta_x(2) .lt. atomic_tol(2)) .and. &
                    (delta_x(3) .lt. atomic_tol(3)) &
                ) then
                    found_match = .true.
                    ! if(debug) write(6,'(a)') "Found matching atoms"
                    exit
                endif
            enddo
            if(.not. found_match) then
                exit
            else
                continue
            endif
        enddo
        if(.not. found_match) then
            cycle
        else
            is_symmetric(ic) = .true.
            if(.not. W_match(W_index(ic))) W_match(W_index(ic)) = .true.
            ! if(debug) write(6,'(a)') "Found matching symmetry operator"
        endif
    enddo

    counter = 0
    do iw = 1, n_W
        if(W_match(iw)) then
            counter = counter + 1
            W(1:3,1:3,counter) = W(1:3,1:3,iw)
        else
            continue
        endif
    enddo
    n_W = count(W_match)

    if(debug) then
        write(6,'(a,i4)') "findTranslationalSymmetry: After finding space group n_W = ", n_W
        do iw = 1, n_W
                write(6,'(a,i4)') "findTranslationalSymmetry: W index = ", iw
            do ix = 1, 3
                write(6,'(3i4)') W(ix,1:3,iw)
            enddo
        enddo
    endif

    n_symm_op = count(is_symmetric)
    allocate(symm_op(3,4,n_symm_op))
    if(debug) write(6,'(a,i6)') "findTranslationalSymmetry: Number of symmetry operators =", n_symm_op

    counter = 0
    do ic = 1, n_candidates
        if(.not. is_symmetric(ic)) cycle
        counter = counter + 1
        if(debug) write(6,'(a,2i4)') &
            "findTranslationalSymmetry: W index, candidate index = ", W_index(ic), ic
        do ix = 1, 3
            symm_op(ix,1:3,counter) = W(ix,1:3,W_index(ic))
            symm_op(ix,4,counter) = nint(12.d0*candidate_translational_operator(ix,ic))
            if(symm_op(ix,4,counter) .eq. 12) symm_op(ix,4,counter) = 0
        enddo
    enddo

    write(6,'(a,i6)') "findTranslationalSymmetry: n_symm_op =", n_symm_op
    write(6,'(a)') "findTranslationalSymmetry: symm_op"
    do ic = 1, n_symm_op
        if(debug) write(6,'(a,i6)') "symmetry operator index", ic
        if(debug) write(6,'(4i4)') symm_op(1,1:4,ic)
        if(debug) write(6,'(4i4)') symm_op(2,1:4,ic)
        if(debug) write(6,'(4i4)') symm_op(3,1:4,ic)
    enddo

    deallocate(n_atom_per_species)
    deallocate(candidate_translational_operator, is_symmetric)
    deallocate(W_index)

end subroutine findTranslationalSymmetry

end module findTranslationalSymmetry_m