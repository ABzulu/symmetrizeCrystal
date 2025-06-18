module symmetrizeConventionalCell_m

    implicit none

    public :: symmetrizeConventionalCell

    private

contains

subroutine symmetrizeConventionalCell( &
    lattice_vectors, n_atom, atomic_coordinates, atomic_species_index, &
    n_symm_op, symm_op, n_site_letters, wyckoff_sites, tol, debug &
)
    integer, intent(in)          :: &
        n_atom, atomic_species_index(n_atom), &
        n_symm_op, symm_op(3,4,192), n_site_letters, wyckoff_sites(3,4,192)
    double precision, intent(in) :: &
        lattice_vectors(3,3), atomic_coordinates(3,n_atom), tol
    logical, intent(in)          :: debug

    integer              :: ia, ja, io, ix, n_sites_found, is
    integer, allocatable :: &
        site_index(:), n_atoms_per_site(:), atom_index_in_site(:,:), counters(:)
    double precision     :: image(3), distance
    logical, allocatable :: assigned(:)

    allocate(site_index(n_atom), assigned(n_atom))

    ! Check which atoms are grouped into same sites
    do io = 1, n_symm_op
        do ia = 1, n_atom
            do ix = 1, 3
                image(ix) = &
                    dble(symm_op(ix,1,io)) * atomic_coordinates(1,ia) + &
                    dble(symm_op(ix,2,io)) * atomic_coordinates(2,ia) + &
                    dble(symm_op(ix,3,io)) * atomic_coordinates(3,ia) + &
                    dble(symm_op(ix,4,io)) / 12.d0
                image(ix) = modulo(image(ix)+0.5d0,1.d0) - 0.5d0
            enddo            

            do ja = 1, n_atom
                if(atomic_species_index(ia) .ne. atomic_species_index(ja)) cycle

                distance = 0
                do ix = 1, 3
                    distance = &
                        distance + &
                        (image(ix) - atomic_coordinates(ix,ja)) * &
                        (image(ix) - atomic_coordinates(ix,ja))
                enddo

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
                        write(6,'(a)') "symmetrizeConventionalCell: Site index mismatch; Stopping program"
                        stop
                    endif
                endif

            enddo
        enddo
    enddo

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
        write(6,'(a,i4)') "symmetrizeConventionalCell: Number of sites found = ", n_sites_found
        do is = 1, n_sites_found
            write(6,'(a,i4)') "symmetrizeConventionalCell: Site index = ", is
            write(6,'(a)') "symmetrizeConventionalCell: Atoms in site = " 
            write(6,'(5i4)') atom_index_in_site(:,is)
        enddo
    endif

    deallocate(site_index, assigned, n_atoms_per_site, atom_index_in_site, counters)

end subroutine symmetrizeConventionalCell

end module symmetrizeConventionalCell_m