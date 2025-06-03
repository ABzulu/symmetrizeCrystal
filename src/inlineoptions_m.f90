module inlineoptions_m
    use m_getopts

    implicit none

    public :: getInlineOptions

    private

contains

subroutine getInlineOptions(output_filename, debug)
    character(len=132), intent(out) :: output_filename
    logical, intent(out) :: debug 

    integer :: nopts, iostat, nargs, nlabels, exst
    character(len=132) :: optArg
    character(len=8) :: optName
    logical, external :: makedirqq

    output_filename = "symmetrizedStructure.out"
    debug = .false.

    nopts = 0
    do
        call getopts('hdo:', optName, optArg, nopts, iostat)
        if( iostat /= 0 ) exit
        select case( optName )
        case( 'h' )
            call manual()
            stop
        case( 'd' )
            debug = .true.
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
    print *, '  -o            Symmetrized output filename'
    print *, '----------------------------------------------------------------'
end subroutine manual

end module inlineoptions_m