module input_m
    use constants, only: ANGSTROM_AU

    implicit none

    public :: readInput, formatInput

    private

contains

subroutine readInput( &
        input_filename, lattice_constant, lattice_vector, &
        n_atom, atomic_coordinates_format, atomic_coordinates, &
        n_species, atomic_species_index &
)
    use fdf, only: &
        fdf_init, fdf_block, fdf_double, fdf_integer, fdf_physical, fdf_string

    character(132), intent(in) :: input_filename
    character(132), intent(out) :: atomic_coordinates_format 
    integer, intent(out) :: n_atom, n_species
    integer, allocatable, intent(out) :: atomic_species_index(:)
    double precision, intent(out) :: lattice_constant, lattice_vector(3,3)
    double precision, allocatable, intent(out) :: atomic_coordinates(:,:)

    integer :: i, ia, iounit

    call fdf_init(input_filename, 'input.fdf.tmp')

    ! This automatically changes the units of lattice constant to Bohr
    lattice_constant = fdf_physical('LatticeConstant',0.d0,'Bohr')
    if(fdf_block('LatticeVectors', iounit)) then
        do i = 1, 3
            read(iounit,*) lattice_vector(i,1:3)
        Enddo
    endif

    atomic_coordinates_format = fdf_string('AtomicCoordinatesFormat','Bohr')
    n_atom = fdf_integer('NumberOfAtoms',0)
    if(n_atom .eq. 0) then
        write(6,*) "Number of atoms = 0; Stopping program"
        stop
    endif
    n_species = fdf_integer('NumberOfSpecies',0)
    if(n_species .eq. 0) then
        write(6,*) "Number of species = 0; Stopping program"
        stop
    endif

    allocate(atomic_coordinates(3,n_atom), atomic_species_index(n_atom))
    if(fdf_block('AtomicCoordinatesAndAtomicSpecies',iounit)) then
        do ia = 1, n_atom
            read(iounit,*) atomic_coordinates(1:3,ia), atomic_species_index(ia)
        enddo
    endif

    do ia = 1, n_atom
        write(6,'(3f16.9,i4)') atomic_coordinates(1:3,ia), atomic_species_index(ia)
    enddo 

end subroutine readInput

subroutine formatInput( &
    lattice_constant, lattice_vector, &
    atomic_coordinates_format, atomic_coordinates, n_atom, debug &
)
    character(132), intent(in) :: atomic_coordinates_format 
    integer, intent(in) :: n_atom
    double precision, intent(in) :: lattice_constant
    double precision, intent(inout) :: lattice_vector(3,3)
    double precision, intent(inout) :: atomic_coordinates(3,n_atom)
    logical, intent(in) :: debug

    integer :: i, ia
    logical :: leqi
    double precision :: atomic_coordinate_temp(3)

    lattice_vector = lattice_vector * lattice_constant

    if(leqi(atomic_coordinates_format,'NotScaledCartesianBohr') .or. &
       leqi(atomic_coordinates_format,'Bohr')) then
        continue
    elseif(leqi(atomic_coordinates_format,'NotScaledCartesianAng') .or. &
           leqi(atomic_coordinates_format,'Ang')) then
        atomic_coordinates = atomic_coordinates * ANGSTROM_AU
    elseif(leqi(atomic_coordinates_format,'ScaledCartesian')) then
        if(lattice_constant .eq. 0.d0) then
            write(6,*) "Lattice constant not given when atomic coordinates format is 'ScaledCartesian'; Stoping program"
            stop
        endif
        atomic_coordinates = atomic_coordinates * lattice_constant
    elseif(leqi(atomic_coordinates_format,'ScaledByLatticeVectors') .or. &
           leqi(atomic_coordinates_format,'Fractional')) then
        do ia = 1, n_atom
            atomic_coordinate_temp(1:3) = atomic_coordinates(1:3,ia)
            do i = 1, 3
                atomic_coordinates(i,ia) = &
                    lattice_vector(i,1) * atomic_coordinate_temp(1) + &
                    lattice_vector(i,2) * atomic_coordinate_temp(2) + &
                    lattice_vector(i,3) * atomic_coordinate_temp(3)
            enddo
        enddo
    else
        write(6,*) "No AtomicCoordinatesFormat given; Stopping program"
        stop
    endif

    if(debug) then
        write(6,*) "Atomic coordinates after formatting"
        do ia = 1, n_atom
            write(6,*) atomic_coordinates(1:3,ia)
        enddo
    endif

end subroutine formatInput

end module input_m