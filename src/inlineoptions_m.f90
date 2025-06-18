module inlineoptions_m
    use m_getopts

    implicit none

    public :: getInlineOptions

    private

contains

subroutine getInlineOptions( &
    lattice_tol, atomic_tol, input_filename, output_filename, maxIteration, debug &
)
    integer, intent(out) :: maxIteration
    double precision, intent(out) :: lattice_tol, atomic_tol
    character(len=132), intent(out) :: input_filename, output_filename
    logical, intent(out) :: debug 

    integer :: nopts, iostat, nargs, nlabels
    character(len=132) :: optArg
    character(len=8) :: optName
    logical, external :: makedirqq

    maxIteration = 100
    lattice_tol = 1.0d-2
    atomic_tol = 1.0d-4
    output_filename = "symmetrizedStructure.out"
    debug = .false.

    nopts = 0
    do
        call getopts('hdi:l:a:o:', optName, optArg, nopts, iostat)
        if( iostat /= 0 ) exit
        select case( optName )
        case( 'h' )
            call manual()
            stop
        case( 'd' )
            debug = .true.
        case( 'i' )
            read(optArg,*) maxIteration
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
    if( nlabels /= 1) then
        print *, 'Use -h option for manual'
        print *, ''
        call manual()
        stop
    endif

    call get_command_argument(nopts,value=input_filename,status=iostat)
    input_filename = trim(input_filename)
    if(iostat /= 0) then
      stop 'Cannot get INPUT file'
    endif

end subroutine getInlineOptions

subroutine manual()
    print *, '----------------------------------------------------------------'
    print *, '  Usage: symmetrizeCrystal [OPTIONS] [INPUT]'
    print *, '  '
    print *, '  INPUT         SIESTA input filename'
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