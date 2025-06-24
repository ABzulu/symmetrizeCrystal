module symmetrizeStructure_m
    use constants, only: eps16

    implicit none

    public :: symmetrizeStructure

    private

contains

subroutine symmetrizeStructure( &
    lattice_vectors, n_atom, atomic_coordinates, atomic_species_index, &
    n_symmetry_operators, symmetry_operators, tol, debug &
)
    integer, intent(in)          :: &
        n_atom, atomic_species_index(n_atom), &
        n_symmetry_operators, symmetry_operators(3,4,192)
    double precision, intent(in) :: lattice_vectors(3,3), tol(3)
    double precision, intent(inout) :: atomic_coordinates(3,n_atom)
    logical, intent(in)          :: debug

    integer              :: counter, ix, iy, iz, ia
    integer, allocatable :: atom_index(:)

    integer                       :: &
        iia, ja, jja, io, n_sites_found, is, &
        matching_atom_index, matching_atom_site_index
    integer, allocatable          :: &
        site_index(:), n_atoms_per_site(:), atom_index_in_site(:,:), counters(:), &
        n_matching_images(:)
    double precision              :: &
        distance, temp_distance, lattice_parameter(3), image(3)
    double precision, allocatable :: temp_atomic_coordintes(:,:), images(:,:)
    logical, allocatable          :: assigned(:)

    if(debug) then
        write(6,'(a)') "symmetrizeStructure: Atomic coordinates before symmetrization"
        do ia = 1, n_atom
            write(6,'(3f16.9)') atomic_coordinates(1:3,ia)
        enddo
    endif

    allocate(site_index(n_atom), assigned(n_atom))

    do ix = 1, 3
        lattice_parameter(ix) = &
            sqrt( &
                lattice_vectors(ix,1)*lattice_vectors(ix,1) + &
                lattice_vectors(ix,2)*lattice_vectors(ix,2) + &
                lattice_vectors(ix,3)*lattice_vectors(ix,3) &
            )
    enddo

    site_index(:) = 0
    assigned(:) = .false.
    n_sites_found = 0
    allocate(images(3,n_atom*125), atom_index(n_atom*125))
    
    do iia = 1, n_atom
        do io = 1, n_symmetry_operators
            do ia = 1, n_atom
                atom_index(ia) = ia
                do ix = 1, 3
                    images(ix,ia) = &
                        dble(symmetry_operators(ix,1,io)) * atomic_coordinates(1,ia) + &
                        dble(symmetry_operators(ix,2,io)) * atomic_coordinates(2,ia) + &
                        dble(symmetry_operators(ix,3,io)) * atomic_coordinates(3,ia) + &
                        dble(symmetry_operators(ix,4,io))/12.d0
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

            do ja = 1, n_atom*125
                if(atomic_species_index(iia) .ne. atomic_species_index(atom_index(ja))) cycle
                if( &
                    (abs(images(1,ja) - atomic_coordinates(1,iia)) .lt. tol(1)) .and. &
                    (abs(images(2,ja) - atomic_coordinates(2,iia)) .lt. tol(2)) .and. &
                    (abs(images(3,ja) - atomic_coordinates(3,iia)) .lt. tol(3)) &
                ) then
                    if (.not. assigned(iia) .and. .not. assigned(atom_index(ja))) then
                        n_sites_found = n_sites_found + 1
                        site_index(iia) = n_sites_found
                        site_index(atom_index(ja)) = n_sites_found
                        assigned(iia) = .true.
                        assigned(atom_index(ja)) = .true.
                    elseif (assigned(iia) .and. .not. assigned(atom_index(ja))) then
                        site_index(atom_index(ja)) = site_index(iia)
                        assigned(atom_index(ja)) = .true.
                    elseif (.not. assigned(iia) .and. assigned(atom_index(ja))) then
                        site_index(iia) = site_index(atom_index(ja))
                        assigned(iia) = .true.
                    elseif (site_index(iia) .ne. site_index(atom_index(ja))) then
                        write(6,'(a)') "symmetrizeStructure: Site index mismatch"
                        write(6,'(a,2i6)') "symmetrizeStructure: iia, ja = ", iia, ja
                        write(6,'(a,3f16.9)') "symmetrizeStructure: image = ", images(1:3,counter)
                        write(6,'(a,6i6)') "symmetrizeStructure: site_index(:) = ", site_index(:)
                        write(6,'(a,6L6)') "symmetrizeStructure: assigned(:) =   ", assigned(:)
                        write(6,'(a)') "symmetrizeStructure: Stopping program"
                        stop
                    endif
                endif
            enddo
        enddo
    enddo
    
    if(n_sites_found .gt. n_atom) then
        write(6,'(a)') "symmetrizeStructure: Too much sites found; Stopping program"
        stop
    endif

    allocate(n_atoms_per_site(n_sites_found), counters(n_sites_found))
    do is = 1, n_sites_found
        n_atoms_per_site(is) = count(site_index(1:n_atom) .eq. is)
    enddo
    allocate(atom_index_in_site(maxval(n_atoms_per_site(1:n_sites_found)),n_sites_found))
    atom_index_in_site(:,:) = 0
    counters(:) = 0
    do ia = 1, n_atom
        counters(site_index(ia)) = counters(site_index(ia)) + 1
        atom_index_in_site(counters(site_index(ia)),site_index(ia)) = ia
    enddo

    if(debug) then
        do is = 1, n_sites_found
            write(6,'(a,i4)') "symmetrizeStructure: Site index = ", is
            write(6,'(a)') "symmetrizeStructure: Atoms in site = " 
            write(6,'(5i4)') atom_index_in_site(:,is)
        enddo
    endif

    !!! This is way too scatch. Revise later
    allocate(n_matching_images(maxval(n_atoms_per_site(1:n_sites_found))))
    allocate(temp_atomic_coordintes(3,n_atom))
    temp_atomic_coordintes(:,:) = 0.d0
    do is = 1, n_sites_found
        n_matching_images(:) = 0
        do ia = 1, n_atoms_per_site(is)
            iia = atom_index_in_site(ia,is)
            do io = 1, n_symmetry_operators
                do ix = 1, 3
                    image(ix) = &
                        dble(symmetry_operators(ix,1,io)) * atomic_coordinates(1,iia) + &
                        dble(symmetry_operators(ix,2,io)) * atomic_coordinates(2,iia) + &
                        dble(symmetry_operators(ix,3,io)) * atomic_coordinates(3,iia) + &
                        dble(symmetry_operators(ix,4,io)) / 12.d0
                    image(ix) = modulo(image(ix)+0.5d0,1.d0) - 0.5d0
                    if(image(ix) + 0.5d0 .lt. eps16) image(ix) = 0.5d0
                enddo
                
                distance = huge(1.d0)
                matching_atom_index = 0
                do ja = 1, n_atoms_per_site(is)
                    jja = atom_index_in_site(ja,is)
                    temp_distance = 0.0d0
                    do ix = 1, 3
                        temp_distance = &
                            temp_distance + &
                            (image(ix) - atomic_coordinates(ix,jja)) * &
                            (image(ix) - atomic_coordinates(ix,jja)) * &
                            lattice_parameter(ix) * &
                            lattice_parameter(ix)
                    enddo
                    temp_distance = sqrt(temp_distance)

                    if(temp_distance .lt. distance) then
                        distance = temp_distance
                        matching_atom_index = jja
                        matching_atom_site_index = ja
                    endif
                enddo

                n_matching_images(matching_atom_site_index) = &
                    n_matching_images(matching_atom_site_index) + 1
                temp_atomic_coordintes(1:3,matching_atom_index) = &
                    temp_atomic_coordintes(1:3,matching_atom_index) + image(1:3)
            enddo
        enddo ! All images in a given site set is made and matched here

        do ia = 1, n_atoms_per_site(is)
            iia = atom_index_in_site(ia,is)
            temp_atomic_coordintes(1:3,iia) = &
                temp_atomic_coordintes(1:3,iia) / dble(n_matching_images(ia))
            if(n_matching_images(ia) .eq. 0) then
                write(6,'(a)') "symmetrizeStructure: tried to divide by 0"
            endif   
        enddo
    enddo

    atomic_coordinates(1:3,1:n_atom) = temp_atomic_coordintes(1:3,1:n_atom)

    if(debug) then
        write(6,'(a)') "symmetrizeStructure: Atomic coordinates after symmetrization"
        do ia = 1, n_atom
            write(6,'(3f16.9)') atomic_coordinates(1:3,ia)
        enddo
    endif

    deallocate(site_index, assigned, n_atoms_per_site, atom_index_in_site, counters)
    deallocate(n_matching_images, temp_atomic_coordintes)

end subroutine symmetrizeStructure

end module symmetrizeStructure_m