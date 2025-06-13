module groupWyckoffOrbits_m
    implicit none

    public :: groupWyckoffOrbits

    private

contains

subroutine groupWyckoffOrbits( &
    symm_op, n_symm_op, n_atom, atomic_coordinates, wyckoff_orbits, debug &
)
    integer, intent(in) :: symm_op(3,4,48), n_symm_op, n_atom
    double precision, intent(in) :: atomic_coordinates
    logical, intent(in) :: debug
    integer, allocatable, intent(out) :: wyckoff_orbits(:)
   
    integer :: ia, ja, io

    allocate(wyckoff_orbits(n_atom))

    do io = 1, n_symm_op
        do ia = 1, n_atom
            do ja = ia+1, n_atom
                
            enddo
        enddo
    enddo

end subroutine groupWyckoffOrbits

end module groupWyckoffOrbits_m