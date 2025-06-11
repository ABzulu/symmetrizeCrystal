module primitiveToConventional_m
    use constants, only: eps12
    use gcd_m, only: gcd

    implicit none

    public :: primitiveToConventional

    private

contains

subroutine primitiveToConventional( &
    crystal_type, W, n_W, debug &
)
    character(len=132), intent(in) :: crystal_type
    integer, intent(in) :: W(3,3,48), n_W
    logical, intent(in) :: debug

    character(len=132) :: junk, crystal_class
    integer :: &
        in, ix, jx, kx, M(3,3), primary_W_index, S(3,3), copy_W(3,3,48), &
        temp_W(3,3), temp_S(3,3), test_e(3), determinant, correction_matrix(3,3)
    double precision :: obv_M(3,3), rev_M(3,3), dbl_M(3,3), temp_dbl_M(3,3)
    logical :: found_e, is_int

    read(crystal_type,*) junk, crystal_class
    crystal_class = trim(crystal_class)

    ! Depending on the Laue class determine conventional cell
    ! 1. Determine primary, secondary, ternary rotation matrices.
    if((crystal_class .eq. '1') .or. (crystal_class .eq. '1_bar')) then
        !!! Make Niggli cell as conventional cell
        continue
    elseif( &
        (crystal_class .eq. '2') .or. &
        (crystal_class .eq. 'm') .or. &
        (crystal_class .eq. '2/m') &
    ) then

        call findPrimaryAxis(W, n_W, 2, M(2,1:3), primary_W_index, debug)

        temp_W(1,1:3) = [1, 0, 0] 
        temp_W(2,1:3) = [0, 1, 0] 
        temp_W(3,1:3) = [0, 0, 1] 
        do in = 1, 2
            do ix = 1, 3
                do jx = 1, 3
                    temp_S(ix,jx) = &
                        temp_W(ix,1)*W(1,jx,primary_W_index) + &
                        temp_W(ix,2)*W(2,jx,primary_W_index) + &
                        temp_W(ix,3)*W(3,jx,primary_W_index)
                enddo
            enddo
            temp_W(1:3,1:3) = temp_S(1:3,1:3)
            do ix = 1, 3
                do jx = 1, 3
                    S(ix,jx) = S(ix,jx) + temp_S(ix,jx)
                enddo
            enddo
        enddo

        if(debug) then
            write(6,'(a)') "primitiveToConventional: S matrix"
            write(6,'(3i4)') S(1,1:3)
            write(6,'(3i4)') S(2,1:3)
            write(6,'(3i4)') S(3,1:3)
        endif

        ! This is not a full proof method but will probably work.
        found_e = .false.
        do ix = -10,10
            do jx = -10,10
                do kx = -10,10
                    test_e(1) = ix
                    test_e(2) = jx
                    test_e(3) = kx
                    if(.not.((S(1,1)*test_e(1)+S(1,2)*test_e(2)+S(1,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(2,1)*test_e(1)+S(2,2)*test_e(2)+S(2,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(3,1)*test_e(1)+S(3,2)*test_e(2)+S(3,3)*test_e(3)) .eq. 0)) cycle
                    found_e = .true.
                    M(3,1:3) = test_e(1:3)
                    exit
                enddo
                if(found_e) exit
            enddo
            if(found_e) exit
        enddo
        if(.not. found_e) then
            write(6,'(a)') "primitiveToConventional: Could not find appropriate e"
            write(6,'(a)') "primitiveToConventional: Stopping program"
            stop
        endif

        M(3,1:3) = M(3,1:3) / gcd(gcd(abs(M(3,1)),abs(M(3,2))),abs(M(3,3)))
        if(.not. (M(3,1) .eq. 0)) then
            if(M(3,1) .le. 0) M(3,1:3) = -M(3,1:3)
        elseif(.not. (M(3,2) .eq. 0)) then
            if(M(3,2) .le. 0) M(3,1:3) = -M(3,1:3)
        elseif(.not. (M(3,3) .eq. 0)) then
            if(M(3,3) .le. 0) M(3,1:3) = -M(3,1:3)
        else
            write(6,'(a)') "secondary e is 0 vector; Stopping program"
            stop
        endif

        if(debug) then
            write(6,'(a,3i4)') "primitiveToConventional: secondary e = ", M(3,1:3)
        endif

        copy_W(:,:,:) = W(:,:,:)
        ! Make it improper so it skips over
        copy_W(1:3,1:3,primary_W_index) = 0
        do ix = 1, 3
            copy_W(ix,ix,primary_W_index) = -1
        enddo        

        call findPrimaryAxis(W, n_W, 2, M(1,1:3), primary_W_index, debug)

        temp_W(1,1:3) = [1, 0, 0] 
        temp_W(2,1:3) = [0, 1, 0] 
        temp_W(3,1:3) = [0, 0, 1] 
        do in = 1, 2
            do ix = 1, 3
                do jx = 1, 3
                    temp_S(ix,jx) = &
                        temp_W(ix,1)*W(1,jx,primary_W_index) + &
                        temp_W(ix,2)*W(2,jx,primary_W_index) + &
                        temp_W(ix,3)*W(3,jx,primary_W_index)
                enddo
            enddo
            temp_W(1:3,1:3) = temp_S(1:3,1:3)
            do ix = 1, 3
                do jx = 1, 3
                    S(ix,jx) = S(ix,jx) + temp_S(ix,jx)
                enddo
            enddo
        enddo

        if(debug) then
            write(6,'(a)') "primitiveToConventional: S matrix"
            write(6,'(3i4)') S(1,1:3)
            write(6,'(3i4)') S(2,1:3)
            write(6,'(3i4)') S(3,1:3)
        endif

        ! This is not a full proof method but will probably work.
        found_e = .false.
        do ix = -10,10
            do jx = -10,10
                do kx = -10,10
                    test_e(1) = ix
                    test_e(2) = jx
                    test_e(3) = kx
                    if(.not.((S(1,1)*test_e(1)+S(1,2)*test_e(2)+S(1,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(2,1)*test_e(1)+S(2,2)*test_e(2)+S(2,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(3,1)*test_e(1)+S(3,2)*test_e(2)+S(3,3)*test_e(3)) .eq. 0)) cycle
                    found_e = .true.
                    M(1,1:3) = test_e(1:3)
                    exit
                enddo
                if(found_e) exit
            enddo
            if(found_e) exit
        enddo
        if(.not. found_e) then
            write(6,'(a)') "primitiveToConventional: Could not find appropriate e"
            write(6,'(a)') "primitiveToConventional: Stopping program"
            stop
        endif

        M(3,1:3) = M(3,1:3) / gcd(gcd(abs(M(3,1)),abs(M(3,2))),abs(M(3,3)))
        if(.not. (M(3,1) .eq. 0)) then
            if(M(3,1) .le. 0) M(3,1:3) = -M(3,1:3)
        elseif(.not. (M(3,2) .eq. 0)) then
            if(M(3,2) .le. 0) M(3,1:3) = -M(3,1:3)
        elseif(.not. (M(3,3) .eq. 0)) then
            if(M(3,3) .le. 0) M(3,1:3) = -M(3,1:3)
        else
            write(6,'(a)') "secondary e is 0 vector; Stopping program"
            stop
        endif

        if(debug) then
            write(6,'(a,3i4)') "primitiveToConventional: secondary e = ", M(3,1:3)
        endif

        do ix = 1, 3
            M(2,ix) = &
                W(ix,1,primary_W_index)*M(1,1) + &
                W(ix,2,primary_W_index)*M(1,2) + &
                W(ix,3,primary_W_index)*M(1,3)
        enddo

        if((crystal_class .eq. '2/m')) then
            determinant = &
                M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) + &
                M(1,2)*(M(2,3)*M(3,1) - M(2,1)*M(3,3)) + &
                M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
            if( .not. (determinant .eq. 2)) then
                write(6,'(a)') "primitiveToConventional: Determinant of M is not 2 when it should be"        
                write(6,'(a)') "primitiveToConventional: Stopping program"
                stop
            endif
            correction_matrix(1,1:3) = [ 1, 0, 0]
            correction_matrix(2,1:3) = [ 0, 1, 0]
            correction_matrix(3,1:3) = [ 0, 0, 1]  
            do ix = 1, 3
                if(((abs(M(ix,1)) + abs(M(ix,2)) + abs(M(ix,3))) .eq. 1)) then
                    if((abs(M(ix,1)) .eq. 1)) then
                        correction_matrix(1,1:3) = [ 0, 0, 1]
                        correction_matrix(2,1:3) = [ 0,-1, 0]
                        correction_matrix(3,1:3) = [ 1, 0, 0]
                        if(debug) then
                            write(6,'(a)') "primitiveToConventional: Centering type is A"
                        endif
                        exit
                    elseif((abs(M(ix,2)) .eq. 1)) then
                        correction_matrix(1,1:3) = [ 0, 1, 0]
                        correction_matrix(2,1:3) = [ 0, 0, 1]
                        correction_matrix(3,1:3) = [ 1, 0, 0]
                        if(debug) then
                            write(6,'(a)') "primitiveToConventional: Centering type is B"
                        endif
                        exit
                    elseif((abs(M(ix,2)) .eq. 1)) then
                        if(debug) then
                            write(6,'(a)') "primitiveToConventional: Centering type is C"
                        endif
                        exit
                    endif
                endif
            enddo
            if( &
                ((abs(M(1,1)) + abs(M(1,2)) + abs(M(1,3))) .eq. 2) .and. &
                ((abs(M(2,1)) + abs(M(2,2)) + abs(M(2,3))) .eq. 2) .and. &
                ((abs(M(3,1)) + abs(M(3,2)) + abs(M(3,3))) .eq. 2) &
            ) then
                correction_matrix(1,1:3) = [ 1, 0,-1]
                correction_matrix(2,1:3) = [ 0, 1, 0]
                correction_matrix(3,1:3) = [ 1, 0, 0]                
            endif

            temp_W(1:3,1:3) = M(1:3,1:3)
            do ix = 1, 3
                do jx = 1, 3
                    M(ix,jx) = &
                        temp_W(ix,1)*correction_matrix(1,jx) + &
                        temp_W(ix,2)*correction_matrix(2,jx) + &
                        temp_W(ix,3)*correction_matrix(3,jx)
                enddo
            enddo
        endif

        if(debug) then
            write(6,'(a)') "primitiveToConventional: M matrix"
            write(6,'(3i4)') M(1,1:3)
            write(6,'(3i4)') M(2,1:3)
            write(6,'(3i4)') M(3,1:3)
        endif

    elseif( &
        (crystal_class .eq. '4') .or. &
        (crystal_class .eq. '4_bar') .or. &
        (crystal_class .eq. '4/m') .or. &
        (crystal_class .eq. '422') .or. &
        (crystal_class .eq. '4mm') .or. &
        (crystal_class .eq. '4_bar2m') .or. &
        (crystal_class .eq. '4/mmm') &
    ) then

        call findPrimaryAxis(W, n_W, 4, M(3,1:3), primary_W_index, debug)

        temp_W(1,1:3) = [1, 0, 0] 
        temp_W(2,1:3) = [0, 1, 0] 
        temp_W(3,1:3) = [0, 0, 1] 
        do in = 1, 4
            do ix = 1, 3
                do jx = 1, 3
                    temp_S(ix,jx) = &
                        temp_W(ix,1)*W(1,jx,primary_W_index) + &
                        temp_W(ix,2)*W(2,jx,primary_W_index) + &
                        temp_W(ix,3)*W(3,jx,primary_W_index)
                enddo
            enddo
            temp_W(1:3,1:3) = temp_S(1:3,1:3)
            do ix = 1, 3
                do jx = 1, 3
                    S(ix,jx) = S(ix,jx) + temp_S(ix,jx)
                enddo
            enddo
        enddo

        if(debug) then
            write(6,'(a)') "primitiveToConventional: S matrix"
            write(6,'(3i4)') S(1,1:3)
            write(6,'(3i4)') S(2,1:3)
            write(6,'(3i4)') S(3,1:3)
        endif

        ! This is not a full proof method but will probably work.
        found_e = .false.
        do ix = -10,10
            do jx = -10,10
                do kx = -10,10
                    test_e(1) = ix
                    test_e(2) = jx
                    test_e(3) = kx
                    if(.not.((S(1,1)*test_e(1)+S(1,2)*test_e(2)+S(1,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(2,1)*test_e(1)+S(2,2)*test_e(2)+S(2,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(3,1)*test_e(1)+S(3,2)*test_e(2)+S(3,3)*test_e(3)) .eq. 0)) cycle
                    found_e = .true.
                    M(1,1:3) = test_e(1:3)
                    exit
                enddo
                if(found_e) exit
            enddo
            if(found_e) exit
        enddo
        if(.not. found_e) then
            write(6,'(a)') "primitiveToConventional: Could not find appropriate e"
            write(6,'(a)') "primitiveToConventional: Stopping program"
            stop
        endif

        M(1,1:3) = M(1,1:3) / gcd(gcd(abs(M(1,1)),abs(M(1,2))),abs(M(1,3)))
        if(.not. (M(1,1) .eq. 0)) then
            if(M(1,1) .le. 0) M(1,1:3) = -M(1,1:3)
        elseif(.not. (M(1,2) .eq. 0)) then
            if(M(1,2) .le. 0) M(1,1:3) = -M(1,1:3)
        elseif(.not. (M(1,3) .eq. 0)) then
            if(M(1,3) .le. 0) M(1,1:3) = -M(1,1:3)
        else
            write(6,'(a)') "secondary e is 0 vector; Stopping program"
            stop
        endif

        if(debug) then
            write(6,'(a,3i4)') "primitiveToConventional: secondary e = ", M(1,1:3)
        endif

        if(.not. ((M(1,1)*M(3,1) + M(1,2)*M(3,2) + M(1,3)*M(3,3)) .lt. eps12)) then
            write(6,'(a)') "primitiveToConventional: primary and ternary vectors are not orthogonal"
            write(6,'(a)') "primitiveToConventional: At 1 Stopping program"
            stop
        endif
        do ix = 1, 3
            M(2,ix) = &
                W(ix,1,primary_W_index)*M(1,1) + &
                W(ix,2,primary_W_index)*M(1,2) + &
                W(ix,3,primary_W_index)*M(1,3)
        enddo

        if(debug) then
            write(6,'(a)') "primitiveToConventional: M matrix"
            write(6,'(3i4)') M(1,1:3)
            write(6,'(3i4)') M(2,1:3)
            write(6,'(3i4)') M(3,1:3)
        endif

    elseif( &
        (crystal_class .eq. '3') .or. &
        (crystal_class .eq. '3_bar') .or. &
        (crystal_class .eq. '32') .or. &
        (crystal_class .eq. '3m') .or. &
        (crystal_class .eq. '3_barm') .or. &
        (crystal_class .eq. '6') .or. &
        (crystal_class .eq. '6_bar') .or. &
        (crystal_class .eq. '6/m') .or. &
        (crystal_class .eq. '622') .or. &
        (crystal_class .eq. '6mm') .or. &
        (crystal_class .eq. '6_bar2m') .or. &
        (crystal_class .eq. '6/mmm') &
    ) then

        call findPrimaryAxis(W, n_W, 3, M(3,1:3), primary_W_index, debug)

        temp_W(1,1:3) = [1, 0, 0] 
        temp_W(2,1:3) = [0, 1, 0] 
        temp_W(3,1:3) = [0, 0, 1] 
        do in = 1, 3
            do ix = 1, 3
                do jx = 1, 3
                    temp_S(ix,jx) = &
                        temp_W(ix,1)*W(1,jx,primary_W_index) + &
                        temp_W(ix,2)*W(2,jx,primary_W_index) + &
                        temp_W(ix,3)*W(3,jx,primary_W_index)
                enddo
            enddo
            temp_W(1:3,1:3) = temp_S(1:3,1:3)
            do ix = 1, 3
                do jx = 1, 3
                    S(ix,jx) = S(ix,jx) + temp_S(ix,jx)
                enddo
            enddo
        enddo

        if(debug) then
            write(6,'(a)') "primitiveToConventional: S matrix"
            write(6,'(3i4)') S(1,1:3)
            write(6,'(3i4)') S(2,1:3)
            write(6,'(3i4)') S(3,1:3)
        endif

        ! This is not a full proof method but will probably work.
        found_e = .false.
        do ix = -10,10
            do jx = -10,10
                do kx = -10,10
                    test_e(1) = ix
                    test_e(2) = jx
                    test_e(3) = kx
                    if(.not.((S(1,1)*test_e(1)+S(1,2)*test_e(2)+S(1,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(2,1)*test_e(1)+S(2,2)*test_e(2)+S(2,3)*test_e(3)) .eq. 0)) cycle
                    if(.not.((S(3,1)*test_e(1)+S(3,2)*test_e(2)+S(3,3)*test_e(3)) .eq. 0)) cycle
                    found_e = .true.
                    M(1,1:3) = test_e(1:3)
                    exit
                enddo
                if(found_e) exit
            enddo
            if(found_e) exit
        enddo
        if(.not. found_e) then
            write(6,'(a)') "primitiveToConventional: Could not find appropriate e"
            write(6,'(a)') "primitiveToConventional: Stopping program"
            stop
        endif

        M(1,1:3) = M(1,1:3) / gcd(gcd(abs(M(1,1)),abs(M(1,2))),abs(M(1,3)))
        if(.not. (M(1,1) .eq. 0)) then
            if(M(1,1) .le. 0) M(1,1:3) = -M(1,1:3)
        elseif(.not. (M(1,2) .eq. 0)) then
            if(M(1,2) .le. 0) M(1,1:3) = -M(1,1:3)
        elseif(.not. (M(1,3) .eq. 0)) then
            if(M(1,3) .le. 0) M(1,1:3) = -M(1,1:3)
        else
            write(6,'(a)') "secondary e is 0 vector; Stopping program"
            stop
        endif

        if(debug) then
            write(6,'(a,3i4)') "primitiveToConventional: secondary e = ", M(1,1:3)
        endif

        if(.not. ((M(1,1)*M(3,1) + M(1,2)*M(3,2) + M(1,3)*M(3,3)) .lt. eps12)) then
            write(6,'(a)') "primitiveToConventional: primary and ternary vectors are not orthogonal"
            write(6,'(a)') "primitiveToConventional: At 1 Stopping program"
            stop
        endif
        do ix = 1, 3
            M(2,ix) = &
                W(ix,1,primary_W_index)*M(1,1) + &
                W(ix,2,primary_W_index)*M(1,2) + &
                W(ix,3,primary_W_index)*M(1,3)
        enddo

        if(.not. ((M(1,1)*M(3,1) + M(1,2)*M(3,2) + M(1,3)*M(3,3)) .lt. eps12)) then
            write(6,'(a)') "primitiveToConventional: primary and ternary vectors are not orthogonal"
            write(6,'(a)') "primitiveToConventional: At 1 Stopping program"
            stop
        endif

        do ix = 1, 3
            M(2,ix) = &
                W(ix,1,primary_W_index)*M(1,1) + &
                W(ix,2,primary_W_index)*M(1,2) + &
                W(ix,3,primary_W_index)*M(1,3)
        enddo

        if(debug) then
            write(6,'(a)') "primitiveToConventional: M matrix"
            write(6,'(3i4)') M(1,1:3)
            write(6,'(3i4)') M(2,1:3)
            write(6,'(3i4)') M(3,1:3)
        endif

    elseif( &
        (crystal_class .eq. '222') .or. &
        (crystal_class .eq. 'mm2') .or. &
        (crystal_class .eq. 'mmm') .or. &
        (crystal_class .eq. '23') .or. &
        (crystal_class .eq. 'm3_bar') &
    ) then

        call findPrimaryAxis(W, n_W, 2, M(1,1:3), primary_W_index, debug)

        copy_W(:,:,:) = W(:,:,:)

        ! Make it improper so it skips over
        copy_W(1:3,1:3,primary_W_index) = 0
        do ix = 1, 3
            copy_W(ix,ix,primary_W_index) = -1
        enddo

        call findPrimaryAxis(copy_W, n_W, 2, M(2,1:3), primary_W_index, debug)

        ! Make it improper so it skips over
        copy_W(1:3,1:3,primary_W_index) = 0
        do ix = 1, 3
            copy_W(ix,ix,primary_W_index) = -1
        enddo

        call findPrimaryAxis(copy_W, n_W, 2, M(3,1:3), primary_W_index, debug)

        determinant = &
            M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) + &
            M(1,2)*(M(2,3)*M(3,1) - M(2,1)*M(3,3)) + &
            M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
        if(determinant .le. 0) then
            temp_W(1:3,1:3) = M(1:3,1:3)
            M(2,1:3) = M(3,1:3)
            M(3,1:3) = temp_W(2,1:3)
        endif

        if((crystal_class .eq. 'mmm')) then
            if( .not. (determinant .eq. 2)) then
                write(6,'(a)') "primitiveToConventional: Determinant of M is not 2 when it should be"        
                write(6,'(a)') "primitiveToConventional: Stopping program"
                stop
            endif
            correction_matrix(1,1:3) = [ 1, 0, 0]
            correction_matrix(2,1:3) = [ 0, 1, 0]
            correction_matrix(3,1:3) = [ 0, 0, 1]  
            do ix = 1, 3
                if(((abs(M(ix,1)) + abs(M(ix,2)) + abs(M(ix,3))) .eq. 1)) then
                    if((abs(M(ix,1)) .eq. 1)) then
                        correction_matrix(1,1:3) = [ 0, 0, 1]
                        correction_matrix(2,1:3) = [ 0,-1, 0]
                        correction_matrix(3,1:3) = [ 1, 0, 0]
                        if(debug) then
                            write(6,'(a)') "primitiveToConventional: Centering type is A"
                        endif
                        exit
                    elseif((abs(M(ix,2)) .eq. 1)) then
                        correction_matrix(1,1:3) = [ 0, 1, 0]
                        correction_matrix(2,1:3) = [ 0, 0, 1]
                        correction_matrix(3,1:3) = [ 1, 0, 0]
                        if(debug) then
                            write(6,'(a)') "primitiveToConventional: Centering type is B"
                        endif
                        exit
                    elseif((abs(M(ix,2)) .eq. 1)) then
                        if(debug) then
                            write(6,'(a)') "primitiveToConventional: Centering type is C"
                        endif
                        exit
                    endif
                endif
            enddo

            temp_W(1:3,1:3) = M(1:3,1:3)
            do ix = 1, 3
                do jx = 1, 3
                    M(ix,jx) = &
                        temp_W(ix,1)*correction_matrix(1,jx) + &
                        temp_W(ix,2)*correction_matrix(2,jx) + &
                        temp_W(ix,3)*correction_matrix(3,jx)
                enddo
            enddo
        endif        

        if(debug) then
            write(6,'(a)') "primitiveToConventional: M matrix"
            write(6,'(3i4)') M(1,1:3)
            write(6,'(3i4)') M(2,1:3)
            write(6,'(3i4)') M(3,1:3)
        endif

    elseif( &
        (crystal_class .eq. '432') .or. &
        (crystal_class .eq. '4_bar3m') .or. &
        (crystal_class .eq. 'm3_barm') &
    ) then

        call findPrimaryAxis(W, n_W, 4, M(1,1:3), primary_W_index, debug)

        copy_W(:,:,:) = W(:,:,:)
        
        ! Make it improper so it skips over
        copy_W(1:3,1:3,primary_W_index) = 0
        do ix = 1, 3
            copy_W(ix,ix,primary_W_index) = -1
        enddo

        call findPrimaryAxis(copy_W, n_W, 4, M(2,1:3), primary_W_index, debug)

        ! Make it improper so it skips over
        copy_W(1:3,1:3,primary_W_index) = 0
        do ix = 1, 3
            copy_W(ix,ix,primary_W_index) = -1
        enddo

        call findPrimaryAxis(copy_W, n_W, 4, M(3,1:3), primary_W_index, debug)

        determinant = &
            M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) + &
            M(1,2)*(M(2,3)*M(3,1) - M(2,1)*M(3,3)) + &
            M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
        if(determinant .le. 0) then
            temp_W(1:3,1:3) = M(1:3,1:3)
            M(2,1:3) = M(3,1:3)
            M(3,1:3) = temp_W(2,1:3)
        endif

        if(debug) then
            write(6,'(a)') "primitiveToConventional: M matrix"
            write(6,'(3i4)') M(1,1:3)
            write(6,'(3i4)') M(2,1:3)
            write(6,'(3i4)') M(3,1:3)
        endif

    else
        write(6,'(a)') "primitiveToConventional: Wrong crystal type given; Stopping program"
        stop        
    endif

    determinant = &
        M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) + &
        M(1,2)*(M(2,3)*M(3,1) - M(2,1)*M(3,3)) + &
        M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
    ! Rhombohedral cell
    if(determinant .eq. 3) then
        obv_M(1,1:3) = [ 2.d0/3.d0, -1.d0/3.d0, -1.d0/3.d0]
        obv_M(2,1:3) = [ 1.d0/3.d0,  1.d0/3.d0, -2.d0/3.d0]
        obv_M(3,1:3) = [ 1.d0/3.d0,  1.d0/3.d0,  1.d0/3.d0]
        rev_M(1,1:3) = [ 1.d0/3.d0, -2.d0/3.d0,  1.d0/3.d0]
        rev_M(2,1:3) = [ 2.d0/3.d0, -1.d0/3.d0, -1.d0/3.d0]
        rev_M(3,1:3) = [ 1.d0/3.d0,  1.d0/3.d0,  1.d0/3.d0]

        dbl_M(1,1:3) = [ dble(M(1,1)), dble(M(1,2)), dble(M(1,3))]
        dbl_M(2,1:3) = [ dble(M(2,1)), dble(M(2,2)), dble(M(2,3))]
        dbl_M(3,1:3) = [ dble(M(3,1)), dble(M(3,2)), dble(M(3,3))]

        temp_dbl_M(:,:) = dbl_M(:,:)
        do ix = 1, 3
            do jx = 1, 3
                dbl_M(ix,jx) = &
                    temp_dbl_M(ix,1) * obv_M(1,jx) + &
                    temp_dbl_M(ix,2) * obv_M(2,jx) + &
                    temp_dbl_M(ix,3) * obv_M(3,jx)
            enddo
        enddo
        is_int = .true.
        do ix = 1, 3
            do jx = 1, 3
                if(abs(dbl_M(ix,jx) - dble(nint(dbl_M(ix,jx)))) .gt. eps12) then
                    is_int = .false.
                    exit
                endif
            enddo
            if(.not. is_int) exit
        enddo

        if(is_int) then
            M(1,1:3) = [ nint(dbl_M(1,1)), nint(dbl_M(1,2)), nint(dbl_M(1,3))]
            M(2,1:3) = [ nint(dbl_M(2,1)), nint(dbl_M(2,2)), nint(dbl_M(2,3))]
            M(3,1:3) = [ nint(dbl_M(3,1)), nint(dbl_M(3,2)), nint(dbl_M(3,3))]

            if(debug) &
                write(6,'(a)') "primitiveToConventional: Crystal is obv hexagonal"
        else
            dbl_M(1,1:3) = [ dble(M(1,1)), dble(M(1,2)), dble(M(1,3))]
            dbl_M(2,1:3) = [ dble(M(2,1)), dble(M(2,2)), dble(M(2,3))]
            dbl_M(3,1:3) = [ dble(M(3,1)), dble(M(3,2)), dble(M(3,3))]    
            temp_dbl_M(:,:) = dbl_M(:,:)
            do ix = 1, 3
                do jx = 1, 3
                    dbl_M(ix,jx) = &
                        temp_dbl_M(ix,1) * rev_M(1,jx) + &
                        temp_dbl_M(ix,2) * rev_M(2,jx) + &
                        temp_dbl_M(ix,3) * rev_M(3,jx)
                enddo
            enddo
            is_int = .true.
            do ix = 1, 3
                do jx = 1, 3
                    if(abs(dbl_M(ix,jx) - dble(nint(dbl_M(ix,jx)))) .gt. eps12) then
                        is_int = .false.
                        exit
                    endif
                enddo
                if(.not. is_int) exit
            enddo
            
            if(is_int) then
                M(1,1:3) = [ nint(dbl_M(1,1)), nint(dbl_M(1,2)), nint(dbl_M(1,3))]
                M(2,1:3) = [ nint(dbl_M(2,1)), nint(dbl_M(2,2)), nint(dbl_M(2,3))]
                M(3,1:3) = [ nint(dbl_M(3,1)), nint(dbl_M(3,2)), nint(dbl_M(3,3))]

                if(debug) &
                    write(6,'(a)') "primitiveToConventional: Crystal is rev hexagonal"
            else
                write(6,'(a)') "primitiveToConventional: Crystal should not be rhombohedral"
            endif
        endif
    endif

end subroutine primitiveToConventional

subroutine findPrimaryAxis(W, n_W, order, primary_axis, primary_W_index, debug)
    integer, intent(in) :: W(3,3,48), n_W, order
    logical, intent(in) :: debug
    integer, intent(out) :: primary_axis(3), primary_W_index

    integer :: iw, in, ix, jx, determinant, test_W(3,3), temp_W(3,3), temp_M(3)
    logical :: W_is_primary(48)

    W_is_primary(:) = .false.
    do iw = 1, n_W
        determinant = &
            W(1,1,iw)*(W(2,2,iw)*W(3,3,iw) - W(2,3,iw)*W(3,2,iw)) + &
            W(1,2,iw)*(W(2,3,iw)*W(3,1,iw) - W(2,1,iw)*W(3,3,iw)) + &
            W(1,3,iw)*(W(2,1,iw)*W(3,2,iw) - W(2,2,iw)*W(3,1,iw))
        if(determinant .eq. -1) cycle
        test_W(1:3,1:3) = W(1:3,1:3,iw)
        if( &
            (test_W(1,1) .eq. 1) .and. (test_W(1,2) .eq. 0) .and. (test_W(1,3) .eq. 0) .and. &
            (test_W(2,1) .eq. 0) .and. (test_W(2,2) .eq. 1) .and. (test_W(2,3) .eq. 0) .and. &
            (test_W(3,1) .eq. 0) .and. (test_W(3,2) .eq. 0) .and. (test_W(3,3) .eq. 1) &
        ) cycle
        do in = 1, order-1
            temp_W = test_W
            do ix = 1, 3
                do jx = 1, 3
                    test_W(ix,jx) = &
                        temp_W(ix,1) * W(1,jx,iw) + &
                        temp_W(ix,2) * W(2,jx,iw) + &
                        temp_W(ix,3) * W(3,jx,iw)
                enddo
            enddo
            if( &
                (test_W(1,1) .eq. 1) .and. (test_W(1,2) .eq. 0) .and. (test_W(1,3) .eq. 0) .and. &
                (test_W(2,1) .eq. 0) .and. (test_W(2,2) .eq. 1) .and. (test_W(2,3) .eq. 0) .and. &
                (test_W(3,1) .eq. 0) .and. (test_W(3,2) .eq. 0) .and. (test_W(3,3) .eq. 1) &
            ) then
                if(in .eq. order-1) W_is_primary(iw) = .true.
                exit
            endif
        enddo
    enddo

    if(debug) then
        write(6,'(a,i4)') "findPrimaryAxis: Primary matrices with order", order
        do iw = 1, n_W
            if(.not. W_is_primary(iw)) cycle
            write(6,'(a,i4)') "findPrimaryAxis: Matrice index = ", iw
            do ix = 1, 3
                write(6,'(3i4)') W(ix,1:3,iw)
            enddo
        enddo
    endif

    temp_M(:) = 0
    do iw = 1, n_W
        if(.not. W_is_primary(iw)) cycle
        ! Check all W with desired order.
        ! As long as x-component of either one is not zero, compare the two x-components
        ! and choose W with the larger x-component. If both are the same choose one with
        ! positive x-component.
        ! If both is zero move to y, and then z
        if(count(W_is_primary) .gt. 1) then
            temp_M(1:3) = [ &
                W(3,2,iw) - W(2,3,iw), &
                W(1,3,iw) - W(3,1,iw), &
                W(2,1,iw) - W(1,2,iw) ]
            temp_M(1:3) = temp_M(1:3) / gcd(gcd(abs(temp_M(1)),abs(temp_M(2))),abs(temp_M(3)))
            if((.not.(temp_M(1) .eq. 0)) .or. (.not.(primary_axis(1) .eq. 0))) then
                if(abs(temp_M(1)) .gt. abs(primary_axis(1))) then
                    primary_axis(1:3) = temp_M(1:3)
                    primary_W_index = iw
                elseif(abs(temp_M(1)) .eq. abs(primary_axis(1))) then
                    if(temp_M(1) .gt. 0) then
                        primary_axis(1:3) = temp_M(1:3)
                        primary_W_index = iw
                    endif
                endif
            elseif((.not.(temp_M(2) .eq. 0)) .or. (.not.(primary_axis(2) .eq. 0))) then
                if(abs(temp_M(2)) .gt. abs(primary_axis(2))) then
                    primary_axis(1:3) = temp_M(1:3)
                    primary_W_index = iw
                elseif(abs(temp_M(1)) .eq. abs(primary_axis(1))) then
                    if(temp_M(2) .gt. 0) then
                        primary_axis(1:3) = temp_M(1:3)
                        primary_W_index = iw
                    endif
                endif
            elseif((.not.(temp_M(3) .eq. 0)) .or. (.not.(primary_axis(3) .eq. 0))) then
                if(abs(temp_M(3)) .gt. abs(primary_axis(3))) then
                    primary_axis(1:3) = temp_M(1:3)
                    primary_W_index = iw
                elseif(abs(temp_M(1)) .eq. abs(primary_axis(1))) then
                    if(temp_M(3) .gt. 0) then
                        primary_axis(1:3) = temp_M(1:3)
                        primary_W_index = iw
                    endif
                endif
            endif
        else
            primary_W_index = iw
            primary_axis(1:3) = [ &
                W(3,2,iw) - W(2,3,iw), &
                W(1,3,iw) - W(3,1,iw), &
                W(2,1,iw) - W(1,2,iw) ]
            primary_axis(1:3) = &
                primary_axis(1:3) / &
                gcd(gcd(abs(primary_axis(1)),abs(primary_axis(2))),abs(primary_axis(3)))    
        endif
    enddo

    if(debug) &
        write(6,'(a,4i4)') &
            "findPrimaryAxis: Primary Axis, primary W index = ", &
            primary_axis(1:3), primary_W_index

end subroutine findPrimaryAxis

end module primitiveToConventional_m