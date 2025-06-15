module niggliReduction_m
    use constants, only: eps6

    implicit none

    public :: niggliReduction

    private

contains

subroutine niggliReduction( &
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

        if(debug) then
            if(debug) write(6,'(a)') "niggliReduction: Lattice vectors"
            write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
            write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
            write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
        endif

        ! Force a_1 <= a_2 <= a_3
        if(G(2,2) - G(3,3) .gt. lattice_tol) then
            if(debug) write(6,'(a)') "niggliReduction: a_1 <= a_2 <= a_3 # 1"
            temp_array(1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)
        endif
        if(G(1,1) - G(2,2) .gt. lattice_tol) then
            if(debug) write(6,'(a)') "niggliReduction: a_1 <= a_2 <= a_3 # 2"
            temp_array(1:3) = reduced_lattice_vectors(1,1:3)
            reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = temp_array(1:3)
        endif
        if(G(2,2) - G(3,3) .gt. lattice_tol) then
            if(debug) write(6,'(a)') "niggliReduction: a_1 <= a_2 <= a_3 # 3"
            temp_array(1:3) = reduced_lattice_vectors(2,1:3)
            reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
            reduced_lattice_vectors(3,1:3) = temp_array(1:3)
        endif

        ! Tie-breaker
        if( &
            (abs(G(2,2) - G(3,3)) .lt. lattice_tol) .and. &
            (abs(G(1,1) - G(2,2)) .lt. lattice_tol) &
        ) then
            if( &
                (abs(G(2,3) - G(3,1)) .gt. eps6) .and. &
                (abs(G(3,1) - G(1,2)) .gt. eps6) &
            ) then
                if(debug) write(6,'(a)') "niggliReduction: Tie-breaker # 1"
                temp_array(1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
                reduced_lattice_vectors(3,1:3) = temp_array(1:3)
                temp_array(1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(2,1:3) = temp_array(1:3)
            else
                continue                
            endif
        endif
        if( &
            (abs(G(2,2) - G(3,3)) .lt. lattice_tol) .and. &
            (abs(G(1,1) - G(2,2)) .lt. lattice_tol) &
        ) then
            if( &
                (abs(G(2,3) - G(3,1)) .gt. eps6) .and. &
                (abs(G(3,1) - G(1,2)) .gt. eps6) &
            ) then
                if(debug) write(6,'(a)') "niggliReduction: Tie-breaker # 2"
                temp_array(1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
                reduced_lattice_vectors(3,1:3) = temp_array(1:3)
                temp_array(1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(2,1:3) = temp_array(1:3)
            else
                continue                
            endif
        endif
        
        if(abs(G(1,1) - G(2,2)) .lt. lattice_tol) then
            if(abs(G(3,1) - G(1,2)) .gt. eps6) then
                continue
            else
                if(debug) write(6,'(a)') "niggliReduction: Tie-breaker # 3"
                temp_array(1:3) = reduced_lattice_vectors(2,1:3)
                reduced_lattice_vectors(2,1:3) = reduced_lattice_vectors(3,1:3)
                reduced_lattice_vectors(3,1:3) = temp_array(1:3)
            endif
        elseif(abs(G(2,2) - G(3,3)) .lt. lattice_tol) then
            if(abs(G(2,3) - G(3,1)) .gt. eps6) then
                continue
            else
                if(debug) write(6,'(a)') "niggliReduction: Tie-breaker # 4"
                temp_array(1:3) = reduced_lattice_vectors(3,1:3)
                reduced_lattice_vectors(1,1:3) = reduced_lattice_vectors(3,1:3)
                reduced_lattice_vectors(3,1:3) = temp_array(1:3)
            endif
        endif

        ! Force negative or zero inner products between lattice vectors.
        inner_products_are_negetive = .false.
        if((G(1,2) .lt. eps6) .and. (G(2,3) .lt. eps6) .and. (G(3,1) .lt. eps6)) &
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

    if(debug) then
        write(6,'(a)') "niggliReduction: lattice vectors after Niggli reduction"
        write(6,'(3f16.9)') reduced_lattice_vectors(1,1:3)
        write(6,'(3f16.9)') reduced_lattice_vectors(2,1:3)
        write(6,'(3f16.9)') reduced_lattice_vectors(3,1:3)
    endif

end subroutine niggliReduction

end module niggliReduction_m