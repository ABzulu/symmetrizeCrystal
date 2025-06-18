module identifyCrystal_m
    implicit none

    public :: identifyCrystal

    private

contains

subroutine identifyCrystal(W, n_W, W_type, crystal_type, debug)
    integer, intent(in) :: W(3,3,48), n_W
    logical, intent(in) :: debug
    integer, intent(out) :: W_type(48)
    character(len=132), intent(out) :: crystal_type

    integer :: iw, ix, trace, determinant, counter(10)

    write(6,'(a)') "identifyCrystal: Type of rotation matrix W"
    
    counter(:) = 0
    do iw = 1, n_W
        trace = 0
        do ix = 1, 3
            trace = trace + W(ix,ix,iw)
        enddo
        determinant = &
            W(1,1,iw)*(W(2,2,iw)*W(3,3,iw) - W(2,3,iw)*W(3,2,iw)) + &
            W(1,2,iw)*(W(2,3,iw)*W(3,1,iw) - W(2,1,iw)*W(3,3,iw)) + &
            W(1,3,iw)*(W(2,1,iw)*W(3,2,iw) - W(2,2,iw)*W(3,1,iw))
        if(determinant .eq. 1) then
            select case(trace)
            case(3)
                W_type(iw) = 1
                counter(6) = counter(6) + 1
            case(-1)
                W_type(iw) = 2
                counter(7) = counter(7) + 1
            case(0)
                W_type(iw) = 3
                counter(8) = counter(8) + 1
            case(1)
                W_type(iw) = 4
                counter(9) = counter(9) + 1
            case(2)
                W_type(iw) = 6
                counter(10) = counter(10) + 1
            end select
        elseif(determinant .eq. -1) then
            select case(trace)
            case(-3)
                W_type(iw) = -1
                counter(5) = counter(5) + 1
            case(1)
                W_type(iw) = -2
                counter(4) = counter(4) + 1
            case(0)
                W_type(iw) = -3
                counter(3) = counter(3) + 1
            case(-1)
                W_type(iw) = -4
                counter(2) = counter(2) + 1
            case(-2)
                W_type(iw) = -6
                counter(1) = counter(1) + 1
            end select
        else
            write(6,*) "identifyCrystal: W has det that is not 1 or -1; Stopping program"
            stop
        endif

        if(debug) write(6,'(a,2i4)') "identifyCrystal: W index, W type = ", iw, W_type(iw)
    enddo

    write(6,'(a39,10i3)') "identifyCrystal: Type of W            |", -6, -4, -3, -2, -1, 1, 2, 3, 4, 6
    write(6,'(a39,10i3)') "identifyCrystal: Numbers of type of W |", counter(1:10)

    if(.not. (counter(6) .eq. 1)) then
        write(6,'(a)') "identifyCrystal: Number of W with type 1 is not 1; Stopping program"
        stop
    endif

    if( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &
    ) then
        crystal_type = "Triclinic 1"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    )  then
        crystal_type = "Triclinic 1_bar"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Monoclinic 2"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 1) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Monoclinic m"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 1) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Monoclinic 2/m"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Orthorhombic 222"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 2) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Orthorhombic mm2"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 3) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Orthorhombic mmm"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 2) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 4"
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 2) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 4_bar"   
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 2) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 1) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 2) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 4/m" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 5) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 2) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 422" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 4) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 2) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 4mm" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 2) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 2) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 4_bar2m" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 2) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 5) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 5) .and. &
        (counter(8) .eq. 0) .and. (counter(9) .eq. 2) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Tetragonal 4/mmm" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Trigonal 3" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 2) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Trigonal 3_bar" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Trigonal 32" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 3) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Trigonal 3m" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 2) .and. &
        (counter(4) .eq. 3) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Trigonal 3_barm" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 2) &        
    ) then
        crystal_type = "Hexagonal 6" 
    elseif( &
        (counter(1) .eq. 2) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 1) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 0) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Hexagonal 6_bar" 
    elseif( &
        (counter(1) .eq. 2) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 2) .and. &
        (counter(4) .eq. 1) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 2) &        
    ) then
        crystal_type = "Hexagonal 6/m" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 7) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 2) &        
    ) then
        crystal_type = "Hexagonal 622" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 6) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 1) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 2) &        
    ) then
        crystal_type = "Hexagonal 6mm" 
    elseif( &
        (counter(1) .eq. 2) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 4) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Hexagonal 6_bar2m" 
    elseif( &
        (counter(1) .eq. 2) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 2) .and. &
        (counter(4) .eq. 7) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 7) .and. &
        (counter(8) .eq. 2) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 2) &        
    ) then
        crystal_type = "Hexagonal 6/mmm" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 8) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Cubic 23" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 8) .and. &
        (counter(4) .eq. 3) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 8) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Cubic m3_bar" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 0) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 0) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 9) .and. &
        (counter(8) .eq. 8) .and. (counter(9) .eq. 6) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Cubic 432" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 6) .and. (counter(3)  .eq. 0) .and. &
        (counter(4) .eq. 6) .and. (counter(5) .eq. 0) .and. (counter(7)  .eq. 3) .and. &
        (counter(8) .eq. 8) .and. (counter(9) .eq. 0) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Cubic 4_bar3m" 
    elseif( &
        (counter(1) .eq. 0) .and. (counter(2) .eq. 6) .and. (counter(3)  .eq. 8) .and. &
        (counter(4) .eq. 9) .and. (counter(5) .eq. 1) .and. (counter(7)  .eq. 9) .and. &
        (counter(8) .eq. 8) .and. (counter(9) .eq. 6) .and. (counter(10) .eq. 0) &        
    ) then
        crystal_type = "Cubic m3_barm"
    else
        write(6,'(a)') "identifyCrystal: No matching crystal type found; Stopping program"
        stop
    endif

    write(6,'(a,a)') "identifyCrystal: Crystal type = ", trim(crystal_type)

end subroutine identifyCrystal

end module identifyCrystal_m