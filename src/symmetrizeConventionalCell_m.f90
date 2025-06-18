module symmetrizeConventionalCell_m

    implicit none

    public :: symmetrizeConventionalCell

    private

contains

subroutine symmetrizeConventionalCell( &
    n_atom, atomic_coordinates, atomic_species_index, &
    n_symm_op, symm_op, tol, debug &
)
    integer, intent(in)          :: &
        n_atom, atomic_species_index(n_atom), &
        n_symm_op, symm_op(3,4,192)
    double precision, intent(in) :: tol
    double precision, intent(inout) :: atomic_coordinates(3,n_atom)
    logical, intent(in)          :: debug

    integer                       :: &
        ia, iia, ja, jja, io, ix, n_sites_found, is, &
        matching_atom_index, matching_atom_site_index
    integer, allocatable          :: &
        site_index(:), n_atoms_per_site(:), atom_index_in_site(:,:), counters(:), &
        n_matching_images(:)
    double precision              :: image(3), distance, temp_distance
    double precision, allocatable :: temp_atomic_coordintes(:,:)
    logical, allocatable          :: assigned(:)

    if(debug) then
        write(6,'(a)') "symmetrizeConventionalCell: Atomic coordinates"
        do ia = 1, n_atom
            write(6,'(3f16.9)') atomic_coordinates(1:3,ia)
        enddo
    endif

    allocate(site_index(n_atom), assigned(n_atom))

    site_index(:) = 0
    assigned(:) = .false.
    ! Check which atoms are grouped into same sites
    do ia = 1, n_atom
        do ja = 1, n_atom
            if(atomic_species_index(ia) .ne. atomic_species_index(ja)) cycle
            do io = 1, n_symm_op
                do ix = 1, 3
                    image(ix) = &
                        dble(symm_op(ix,1,io)) * atomic_coordinates(1,ia) + &
                        dble(symm_op(ix,2,io)) * atomic_coordinates(2,ia) + &
                        dble(symm_op(ix,3,io)) * atomic_coordinates(3,ia) + &
                        dble(symm_op(ix,4,io)) / 12.d0
                    image(ix) = modulo(image(ix)+0.5d0,1.d0) - 0.5d0
                enddo            

                distance = 0
                do ix = 1, 3
                    distance = &
                        distance + &
                        (image(ix) - atomic_coordinates(ix,ja)) * &
                        (image(ix) - atomic_coordinates(ix,ja))
                enddo
                distance = sqrt(distance)

                if(distance .lt. tol) then
                    if (.not. assigned(ia) .and. .not. assigned(ja)) then
                        n_sites_found = n_sites_found + 1
                        site_index(ia) = n_sites_found
                        site_index(ja) = n_sites_found
                        assigned(ia) = .true.
                        assigned(ja) = .true.
                    elseif (assigned(ia) .and. .not. assigned(ja)) then
                        site_index(ja) = site_index(ia)
                        assigned(ja) = .true.
                    elseif (.not. assigned(ia) .and. assigned(ja)) then
                        site_index(ia) = site_index(ja)
                        assigned(ia) = .true.
                    elseif (site_index(ia) .ne. site_index(ja)) then
                        write(6,'(a)') "symmetrizeConventionalCell: Site index mismatch"
                        write(6,'(a,2i6)') "symmetrizeConventionalCell: ia, ja = ", ia, ja
                        write(6,'(a,3f16.9)') "symmetrizeConventionalCell: image = ", image(1:3)
                        write(6,'(a,6i6)') "symmetrizeConventionalCell: site_index(:) = ", site_index(:)
                        write(6,'(a,6L6)') "symmetrizeConventionalCell: assigned(:) =   ", assigned(:)
                        write(6,'(a)') "symmetrizeConventionalCell: Stopping program"
                        stop
                    endif
                endif

            enddo
        enddo
    enddo

    if(debug) write(6,'(a,i4)') "symmetrizeConventionalCell: Number of sites found = ", n_sites_found
    
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
            write(6,'(a,i4)') "symmetrizeConventionalCell: Site index = ", is
            write(6,'(a)') "symmetrizeConventionalCell: Atoms in site = " 
            write(6,'(5i4)') atom_index_in_site(:,is)
        enddo
    endif

    allocate(n_matching_images(maxval(n_atoms_per_site(1:n_sites_found))))
    allocate(temp_atomic_coordintes(3,n_atom))
    temp_atomic_coordintes(:,:) = 0.d0
    do is = 1, n_sites_found
        n_matching_images(:) = 0
        do ia = 1, n_atoms_per_site(is)
            iia = atom_index_in_site(ia,is)
            do io = 1, n_symm_op
                do ix = 1, 3
                    image(ix) = &
                        dble(symm_op(ix,1,io)) * atomic_coordinates(1,iia) + &
                        dble(symm_op(ix,2,io)) * atomic_coordinates(2,iia) + &
                        dble(symm_op(ix,3,io)) * atomic_coordinates(3,iia) + &
                        dble(symm_op(ix,4,io)) / 12.d0
                    image(ix) = modulo(image(ix)+0.5d0,1.d0) - 0.5d0
                enddo
                
                distance = huge(1.d0)
                matching_atom_index = 0
                do ja = 1, n_atoms_per_site(is)
                    jja = atom_index_in_site(ja,is)
                    temp_distance = 0
                    do ix = 1, 3
                        temp_distance = &
                            temp_distance + &
                            (image(ix) - atomic_coordinates(ix,jja)) * &
                            (image(ix) - atomic_coordinates(ix,jja))
                    enddo
                    temp_distance = sqrt(temp_distance)

                    if(temp_distance .lt. distance) then
                        distance = temp_distance
                        matching_atom_index = jja
                        matching_atom_site_index = ja
                    endif
                enddo

                if(distance .gt. tol) then
                    write(6,'(a)') "symmetrizeConventionalCell: Warning! Distance between image and closest atom larger than tolerance!"
                    write(6,'(a,2f16.9)') "symmetrizeConventionalCell: distance, tol = ", distance, tol
                endif
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
        enddo
    enddo

    atomic_coordinates(1:3,1:n_atom) = temp_atomic_coordintes(1:3,1:n_atom)

    deallocate(site_index, assigned, n_atoms_per_site, atom_index_in_site, counters)
    deallocate(n_matching_images, temp_atomic_coordintes)

end subroutine symmetrizeConventionalCell

end module symmetrizeConventionalCell_m