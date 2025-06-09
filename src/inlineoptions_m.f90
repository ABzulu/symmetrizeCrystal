module inlineoptions_m
    use m_getopts

    implicit none

    public :: getInlineOptions

    private

contains

subroutine getInlineOptions(lattice_tol, atomic_tol, output_filename, debug)
    double precision, intent(out) :: lattice_tol, atomic_tol
    character(len=132), intent(out) :: output_filename
    logical, intent(out) :: debug 

    integer :: nopts, iostat, nargs, nlabels, exst
    character(len=132) :: optArg
    character(len=8) :: optName
    logical, external :: makedirqq

    lattice_tol = 1.0d0
    atomic_tol = 1.0d-5
    output_filename = "symmetrizedStructure.out"
    debug = .false.

    nopts = 0
    do
        call getopts('hdl:a:o:', optName, optArg, nopts, iostat)
        if( iostat /= 0 ) exit
        select case( optName )
        case( 'h' )
            call manual()
            stop
        case( 'd' )
            debug = .true.
        case( 'l' )
            read(optArg,*) lattice_tol
        case( 'a' )
            read(optArg,*) atomic_tol
        case( 'o' )
            read(optArg,*) output_filename
        end select
    enddo
    
    nargs = command_argument_count()
    nlabels = nargs - nopts + 1
    if( nlabels /= 0) then
        print *, 'Use -h option for manual'
        print *, ''
        call manual()
        stop
    endif

end subroutine getInlineOptions

subroutine manual()
    print *, '----------------------------------------------------------------'
    print *, '  Usage: symmetrizeCrystal [OPTIONS]'
    print *, '  '
    print *, '  options: '
    print *, '  '
    print *, '  -h            Print this help'
    print *, '  -d            Enable debug mode'
    print *, '  -l            Lattice tolerance for searching point group'
    print *, '  -a            Atomic tolerance for searching space group'
    print *, '  -o            Symmetrized output filename'
    print *, '----------------------------------------------------------------'
end subroutine manual

end module inlineoptions_m