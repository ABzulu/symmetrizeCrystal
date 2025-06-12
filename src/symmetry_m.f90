module symmetry_m
    use constants, only: eps16

    implicit none

    public :: findRotationalSymmetry, findTranslationalSymmetry

    private

contains

subroutine findRotationalSymmetry(G, lattice_error, W, n_W, debug)
    double precision, intent(in) :: G(3,3), lattice_error
    integer, intent(out) :: W(3,3,48)
    integer, intent(out) :: n_W
    logical, intent(in) :: debug

    integer :: ROT(3,3), det

    integer :: icounter
    double precision :: &
        temp_G(3,3), G_tilda(3,3), diag(3), off_diag(6), cos1, cos2, theta
    integer :: i11, i12, i13, i21, i22, i23, i31, i32, i33, ir, i, j

    W = 0
    n_W = 0
    if(debug) write(6,*) "findRotationalSymmetry: candidate W"
    do i11 = -2,2; do i12 = -2,2; do i13 = -2,2
    do i21 = -2,2; do i22 = -2,2; do i23 = -2,2
    do i31 = -2,2; do i32 = -2,2; do i33 = -2,2
        det = &
            i11*(i22*i33 - i23*i32) - &
            i12*(i21*i33 - i23*i31) + &
            i13*(i21*i32 - i22*i31)
        if((.not.(det .eq. 1)) .and. (.not.(det .eq. -1))) cycle

        ROT(1,1:3) = [i11, i12, i13]
        ROT(2,1:3) = [i21, i22, i23]
        ROT(3,1:3) = [i31, i32, i33]

        ! if(debug) then
        !     write(6,'(3i4)') ROT(1,1:3)
        !     write(6,'(3i4)') ROT(2,1:3)
        !     write(6,'(3i4)') ROT(3,1:3)
        ! endif

        do i = 1, 3
            do j = 1, 3
                temp_G(i,j) = &
                    ROT(1,i) * G(1,j) + &
                    ROT(2,i) * G(2,j) + &
                    ROT(3,i) * G(3,j)
            enddo
        enddo
        do i = 1, 3
            do j = 1, 3
                G_tilda(i,j) = &
                    temp_G(i,1) * ROT(1,j) + &
                    temp_G(i,2) * ROT(2,j) + &
                    temp_G(i,3) * ROT(3,j)
            enddo
        enddo

        ! Diagonal elements represent lengths
        do i = 1, 3
            diag(i) = abs(G_tilda(i,i) - G(i,i))
        enddo
        ! Off-diagonal elements represent angles
        ! https://arxiv.org/html/1808.01590v2 Section V.1 Step(f)
        icounter = 0
        do i = 1, 3
            do j = 1, 3
                if(i == j) cycle
                icounter = icounter + 1
                cos1 = G_tilda(i,j)/sqrt(G_tilda(i,i)*G_tilda(j,j))
                cos2 = G(i,j)/sqrt(G(i,i)*G(j,j))
                theta = &
                    acos(max(-1.d0,min(1.d0,cos1))) - &
                    acos(max(-1.d0,min(1.d0,cos2)))
                off_diag(icounter) = &
                    sin(abs(theta))*sqrt(0.25d0*(G_tilda(i,i)+G(i,i))*(G_tilda(j,j)+G(j,j)))
            enddo
        enddo
        ! if(debug) write(6,*) "findRotationalSymmetry: diag, off_diag"
        ! if(debug) write(6,'(3f16.9)') diag(1:3)
        ! if(debug) write(6,'(3f16.9)') off_diag(1:3)
        ! if(debug) write(6,'(3f16.9)') off_diag(4:6)
        if((maxval(diag) .gt. lattice_error) .or. &
           (maxval(off_diag) .gt. lattice_error)) cycle

        n_W = n_W + 1
        if(n_W .gt. 48) then
            write(6,*) "findRotationalSymmetry: Too many symmetry operations"
            stop
        endif
        W(1:3,1:3,n_W) = ROT(1:3,1:3)

        if(debug) then
            write(6,*) "findRotationalSymmetry: n_W", n_W
            write(6,'(3i6)') W(1,1:3,n_W)
            write(6,'(3i6)') W(2,1:3,n_W)
            write(6,'(3i6)') W(3,1:3,n_W)
        endif
    enddo;enddo;enddo;enddo;enddo;enddo;enddo;enddo;enddo

end subroutine findRotationalSymmetry

subroutine findTranslationalSymmetry( &
    n_atom, atomic_coordinates, n_species, atomic_species_index, &
    atomic_tol, n_W, W, n_symm_op, symm_op, debug &
)
    integer, intent(in) :: n_atom, n_species, atomic_species_index(n_atom)
    double precision, intent(in) :: atomic_coordinates(3,n_atom), atomic_tol
    logical, intent(in) :: debug
    integer, intent(inout) :: n_W, W(3,3,48)
    integer, intent(out) :: n_symm_op
    integer, intent(out) :: symm_op(3,4,192)

    integer :: is, n_candidates, ia, ja, iw, ix, counter, ic
    integer, allocatable :: n_atom_per_species(:), W_index(:)
    double precision :: delta_x(3)
    double precision, allocatable :: candidate_translational_operator(:,:)
    logical :: found_match, W_match(48)
    logical, allocatable :: is_symmetric(:)

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
    if(debug) write(6,*) "findTranslationalSymmetry: candidate_translational_operators"
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
                    if(abs(candidate_translational_operator(ix,counter) - 1.d0) .lt. atomic_tol) &
                        candidate_translational_operator(ix,counter) = 0.d0
                enddo
                if(debug) write(6,'(i6,3f16.9)') &
                    counter, candidate_translational_operator(1:3,counter)
            enddo
        enddo
    enddo

    W_match(:) = .false.
    ! Check if the candidate translational operator maps all the atoms to the
    ! same species of atoms.
    if(debug) write(6,*) "findTranslationalSymmetry: delta_x"
    is_symmetric(:) = .false.
    do ic = 1, n_candidates
        do ia = 1, n_atom
            found_match = .false.
            do ja = 1, n_atom
                if(.not. (atomic_species_index(ia) .eq. atomic_species_index(ja))) cycle                    
                do ix = 1, 3
                    delta_x(ix) = abs( &
                        modulo( &
                            dble(W(ix,1,W_index(ic))) * atomic_coordinates(1,ia) + & 
                            dble(W(ix,2,W_index(ic))) * atomic_coordinates(2,ia) + & 
                            dble(W(ix,3,W_index(ic))) * atomic_coordinates(3,ia) + &
                            candidate_translational_operator(ix,ic) - &
                            atomic_coordinates(ix,ja) + 0.5d0, &
                            1.d0 &
                        ) - 0.5d0&
                    )                       
                enddo
                if(debug) write(6,'(a,3i6,3f16.9)') &
                    "ic, ia, ja, delta_x = ", ic, ia, ja, delta_x(1:3)
                if(maxval(delta_x) .lt. atomic_tol) then
                    found_match = .true.
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
    if(n_symm_op .gt. 192) then
        write(6,'(a)') "findTranslationalSymmetry: Too many symmetry operators found"
        write(6,'(a,i7)') "findTranslationalSymmetry: Number of symmetry operators = ", n_symm_op
        stop
    endif

    counter = 0
    do ic = 1, n_candidates
        if(.not. is_symmetric(ic)) cycle
        counter = counter + 1
        if(debug) write(6,'(a,i4)') "findTranslationalSymmetry: W index = ", W_index(ic)
        do ix = 1, 3
            symm_op(ix,1:3,counter) = W(ix,1:3,W_index(ic))
            symm_op(ix,4,counter) = nint(12.d0*candidate_translational_operator(ix,counter))
        enddo
    enddo

    write(6,*) "findTranslationalSymmetry: n_symm_op =", n_symm_op
    write(6,*) "findTranslationalSymmetry: symm_op"
    do ic = 1, n_symm_op
        write(6,*) "symmetry operator index", ic
        write(6,'(4i4)') symm_op(1,1:4,ic) 
        write(6,'(4i4)') symm_op(2,1:4,ic) 
        write(6,'(4i4)') symm_op(3,1:4,ic) 
    enddo

    deallocate(n_atom_per_species)
    deallocate(candidate_translational_operator, is_symmetric)
    deallocate(W_index)

end subroutine findTranslationalSymmetry

end module symmetry_m