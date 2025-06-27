module cellReduction_m
    use constants, only: eps12

    implicit none

    public :: cellReduction

    private

contains

subroutine cellReduction( &
    lattice_vectors, reduced_lattice_vectors, lattice_tol, debug &
)
    double precision, intent(in) :: lattice_vectors(3,3), lattice_tol
    logical, intent(in) :: debug
    double precision, intent(out) :: reduced_lattice_vectors(3,3)

    integer :: it, ix, jx
    double precision :: G(3,3), temp_array(3)
    logical :: inner_products_are_negetive

    reduced_lattice_vectors(:,:) = lattice_vectors(:,:)

    do it = 1, 100
        do ix = 1, 3;do jx = 1, 3
            G(ix,jx) = &
                reduced_lattice_vectors(ix,1) * reduced_lattice_vectors(jx,1) + &
                reduced_lattice_vectors(ix,2) * reduced_lattice_vectors(jx,2) + &
                reduced_lattice_vectors(ix,3) * reduced_lattice_vectors(jx,3)
        enddo;enddo

        G(1:3,1:3) = G(1:3,1:3) / maxval([G(1,1), G(2,2), G(3,3)])

        if(debug) then
            if(debug) write(6,'(a,i4)') "cellReduction: Lattice vectors for iteration", it
            write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
            write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
            write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
        endif

        if(debug) then
            write(6,'(a)') "cellReduction: Metric G"
            write(6,'(3f16.9)') G(1,1:3)
            write(6,'(3f16.9)') G(2,1:3)
            write(6,'(3f16.9)') G(3,1:3)
        endif

        ! Force a_1 <= a_2 <= a_3
        if((sqrt(G(2,2)) - sqrt(G(3,3))) .gt. lattice_tol) then
            if(debug) write(6,'(a)') "cellReduction: a_1 <= a_2 <= a_3 # 1"
            temp_array(1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)
        endif
        if((sqrt(G(1,1)) - sqrt(G(2,2))) .gt. lattice_tol) then
            if(debug) write(6,'(a)') "cellReduction: a_1 <= a_2 <= a_3 # 2"
            temp_array(1:3) = reduced_lattice_vectors(1,1:3)
            reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = temp_array(1:3)
        endif
        if((sqrt(G(2,2)) - sqrt(G(3,3))) .gt. lattice_tol) then
            if(debug) write(6,'(a)') "cellReduction: a_1 <= a_2 <= a_3 # 3"
            temp_array(1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)
        endif

        ! Force negative or zero inner products between lattice vectors.
        inner_products_are_negetive = .false.
        if((G(1,2) .lt. eps12) .and. (G(2,3) .lt. eps12) .and. (G(3,1) .lt. eps12)) &
            inner_products_are_negetive = .true.

        if(.not. inner_products_are_negetive) then
            if((G(1,2) .gt. G(2,3)) .and. (G(1,2) .gt. G(3,1))) then
                reduced_lattice_vectors(2,1:3) = &
                    reduced_lattice_vectors(2,1:3) - &
                    reduced_lattice_vectors(1,1:3)
            elseif((G(2,3) .gt. G(3,1)) .and. (G(2,3) .gt. G(1,2))) then
                reduced_lattice_vectors(3,1:3) = &
                    reduced_lattice_vectors(3,1:3) - &
                    reduced_lattice_vectors(2,1:3)
            elseif((G(3,1) .gt. G(1,2)) .and. (G(3,1) .gt. G(2,3))) then
                reduced_lattice_vectors(3,1:3) = &
                    reduced_lattice_vectors(3,1:3) - &
                    reduced_lattice_vectors(1,1:3)        
            endif
        endif

        if(inner_products_are_negetive) exit
    enddo

    do ix = 1, 3;do jx = 1, 3
        G(ix,jx) = &
            reduced_lattice_vectors(ix,1) * reduced_lattice_vectors(jx,1) + &
            reduced_lattice_vectors(ix,2) * reduced_lattice_vectors(jx,2) + &
            reduced_lattice_vectors(ix,3) * reduced_lattice_vectors(jx,3)
    enddo;enddo

    G(1:3,1:3) = G(1:3,1:3) / maxval([G(1,1), G(2,2), G(3,3)])

    ! Tie-breaker
    if( &
        abs(sqrt(G(1,1)) - sqrt(G(2,2))) .lt. eps12 .and. &
        abs(sqrt(G(2,2)) - sqrt(G(3,3))) .lt. eps12 .and. &
        abs(sqrt(G(3,3)) - sqrt(G(1,1))) .lt. eps12 &
    ) then
        if( &
            (reduced_lattice_vectors(1,1) .gt. reduced_lattice_vectors(2,1)) .and. &
            (reduced_lattice_vectors(1,1) .gt. reduced_lattice_vectors(3,1)) &
        ) then
            continue
        elseif( &
            (reduced_lattice_vectors(2,1) .gt. reduced_lattice_vectors(3,1)) .and. &
            (reduced_lattice_vectors(2,1) .gt. reduced_lattice_vectors(1,1)) &                
        ) then
            temp_array(1:3) = reduced_lattice_vectors(1,1:3)
            reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = temp_array(1:3)
        elseif( &
            (reduced_lattice_vectors(3,1) .gt. reduced_lattice_vectors(2,1)) .and. &
            (reduced_lattice_vectors(3,1) .gt. reduced_lattice_vectors(1,1)) &                
        ) then
            temp_array(1:3) = reduced_lattice_vectors(1,1:3)
            reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)
        endif
        if(reduced_lattice_vectors(2,2) .gt. reduced_lattice_vectors(3,2)) then
            continue
        else
            temp_array(1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)                
        endif
    elseif(abs(sqrt(G(1,1)) - sqrt(G(2,2))) .lt. eps12) then
        if( &
            (reduced_lattice_vectors(1,1) .gt. reduced_lattice_vectors(2,1)) &
        ) then
            continue
        elseif( &
            (reduced_lattice_vectors(2,1) .gt. reduced_lattice_vectors(1,1)) &
        ) then
            temp_array(1:3) = reduced_lattice_vectors(1,1:3)
            reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = temp_array(1:3)
        endif
    elseif(abs(sqrt(G(2,2)) - sqrt(G(3,3))) .lt. eps12) then
        if( &
            (reduced_lattice_vectors(2,2) .gt. reduced_lattice_vectors(3,2)) &
        ) then
            continue
        elseif( &
            (reduced_lattice_vectors(3,2) .gt. reduced_lattice_vectors(2,2)) &
        ) then
            temp_array(1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)
        endif        
    endif

    if(debug) then
        write(6,'(a)') "cellReduction: Lattice vectors after Niggli reduction"
        write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
        write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
        write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
    endif

end subroutine cellReduction

end module cellReduction_m