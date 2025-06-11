module gcd_m

    implicit none

    public :: gcd

    private

contains

integer function gcd(number1_int,number2_int)
    integer, intent(in) :: number1_int, number2_int
    
    integer :: a, b, temp

    a = number1_int
    b = number2_int
  
    if(number1_int == 0) then
        gcd = number2_int
    elseif(number2_int == 0) then
        gcd = number1_int
    elseif(number1_int < number2_int) then
        temp = a
        a = b
        b = temp    
        do 
            temp = mod(a, b)
            if( temp == 0 ) exit
            a = b
            b = temp
        end do    
        gcd = b
    end if

end function gcd

end module gcd_m