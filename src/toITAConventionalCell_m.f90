module toITAConventionalCell_m
    implicit none

    public :: toITAConventionalCell

    private

contains

subroutine toITAConventionalCell( &
    crystal_type, lattice_vectors, ITA_lattice_vectors, debug &
)
    character(len=132), intent(in) :: crystal_type
    double precision, intent(in) :: lattice_vectors(3,3)
    logical, intent(in) :: debug
    double precision, intent(out) :: ITA_lattice_vectors(3,3)

    character(len=132) :: crystal_system
    integer :: ix
    double precision :: lattice_parameters(3), alpha, beta, gamma

    read(crystal_type,*) crystal_system

    !!! How to distinguish rhombohedral system from hexagonal
    select case(crystal_system)
    case("Triclinic")
        do ix = 1, 3
            lattice_parameters(ix) = &
                lattice_vectors(ix,1) * lattice_vectors(ix,1) + &
                lattice_vectors(ix,2) * lattice_vectors(ix,2) + &
                lattice_vectors(ix,3) * lattice_vectors(ix,3)
            lattice_parameters(ix) = sqrt(lattice_parameters(ix))
        enddo
        alpha = &
            lattice_vectors(2,1) * lattice_vectors(3,1) + &
            lattice_vectors(2,2) * lattice_vectors(3,2) + &
            lattice_vectors(2,3) * lattice_vectors(3,3)
        alpha = &
            alpha / &
            lattice_parameters(2) / &
            lattice_parameters(3)
        beta = &
            lattice_vectors(3,1) * lattice_vectors(1,1) + &
            lattice_vectors(3,2) * lattice_vectors(1,2) + &
            lattice_vectors(3,3) * lattice_vectors(1,3)
        beta = &
            beta / &
            lattice_parameters(3) / &
            lattice_parameters(1)
        gamma = &
            lattice_vectors(1,1) * lattice_vectors(2,1) + &
            lattice_vectors(1,2) * lattice_vectors(2,2) + &
            lattice_vectors(1,3) * lattice_vectors(2,3)
        gamma = &
            gamma / &
            lattice_parameters(1) / &
            lattice_parameters(2)

        ITA_lattice_vectors(1,1:3) = &
            lattice_parameters(1) * [1.d0, 0.d0, 0.d0]
        ITA_lattice_vectors(2,1:3) = &
            lattice_parameters(2) * [cos(beta), sin(beta), 0.d0]
        ITA_lattice_vectors(3,1) = cos(beta)
        ITA_lattice_vectors(3,2) = (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)
        ITA_lattice_vectors(3,3) = &
            sqrt(1 - ITA_lattice_vectors(3,1)*ITA_lattice_vectors(3,1) - &
                ITA_lattice_vectors(3,2)*ITA_lattice_vectors(3,2) &
            )

    case("Monoclinic")
        do ix = 1, 3
            lattice_parameters(ix) = &
                lattice_vectors(ix,1) * lattice_vectors(ix,1) + &
                lattice_vectors(ix,2) * lattice_vectors(ix,2) + &
                lattice_vectors(ix,3) * lattice_vectors(ix,3)
            lattice_parameters(ix) = sqrt(lattice_parameters(ix))
        enddo
        beta = &
            lattice_vectors(3,1) * lattice_vectors(1,1) + &
            lattice_vectors(3,2) * lattice_vectors(1,2) + &
            lattice_vectors(3,3) * lattice_vectors(1,3)
        beta = &
            beta / &
            lattice_parameters(3) / &
            lattice_parameters(1)
        ITA_lattice_vectors(1,1:3) = &
            lattice_parameters(1) * [1.d0, 0.d0, 0.d0]
        ITA_lattice_vectors(2,1:3) = &
            lattice_parameters(2) * [0.d0, 1.d0, 0.d0]
        ITA_lattice_vectors(3,1:3) = &
            lattice_parameters(3) * [0.d0, cos(beta), sin(beta)]
    case("Orthorhombic")
        do ix = 1, 3
            lattice_parameters(ix) = &
                lattice_vectors(ix,1) * lattice_vectors(ix,1) + &
                lattice_vectors(ix,2) * lattice_vectors(ix,2) + &
                lattice_vectors(ix,3) * lattice_vectors(ix,3)
            lattice_parameters(ix) = sqrt(lattice_parameters(ix))
        enddo
        ITA_lattice_vectors(1,1:3) = lattice_parameters(1) * [1.d0, 0.d0, 0.d0]
        ITA_lattice_vectors(2,1:3) = lattice_parameters(2) * [0.d0, 1.d0, 0.d0]                
        ITA_lattice_vectors(3,1:3) = lattice_parameters(3) * [0.d0, 0.d0, 1.d0]                
    case("Tetragonal")
        do ix = 1, 3
            lattice_parameters(ix) = &
                lattice_vectors(ix,1) * lattice_vectors(ix,1) + &
                lattice_vectors(ix,2) * lattice_vectors(ix,2) + &
                lattice_vectors(ix,3) * lattice_vectors(ix,3)
            lattice_parameters(ix) = sqrt(lattice_parameters(ix))
        enddo
        ITA_lattice_vectors(1,1:3) = lattice_parameters(1) * [1.d0, 0.d0, 0.d0]
        ITA_lattice_vectors(2,1:3) = lattice_parameters(2) * [0.d0, 1.d0, 0.d0]                
        ITA_lattice_vectors(3,1:3) = lattice_parameters(3) * [0.d0, 0.d0, 1.d0]
    case("Hexagonal")
        do ix = 1, 3
            lattice_parameters(ix) = &
                lattice_vectors(ix,1) * lattice_vectors(ix,1) + &
                lattice_vectors(ix,2) * lattice_vectors(ix,2) + &
                lattice_vectors(ix,3) * lattice_vectors(ix,3)
            lattice_parameters(ix) = sqrt(lattice_parameters(ix))
        enddo
        ITA_lattice_vectors(1,1:3) = lattice_parameters(1) * [ 1.0d0, 0.0d0, 0.0d0]
        ITA_lattice_vectors(2,1:3) = lattice_parameters(2) * [-0.5d0, sqrt(3.d0)*0.5d0, 0.0d0]                
        ITA_lattice_vectors(3,1:3) = lattice_parameters(3) * [ 0.0d0, 0.0d0, 1.0d0]        
    case("Cubic")
        do ix = 1, 3
            lattice_parameters(ix) = &
                lattice_vectors(1,1) * lattice_vectors(1,1) + &
                lattice_vectors(1,2) * lattice_vectors(1,2) + &
                lattice_vectors(1,3) * lattice_vectors(1,3)
            lattice_parameters(ix) = sqrt(lattice_parameters(ix))
        enddo
        ITA_lattice_vectors(1,1:3) = lattice_parameters(1) * [1.d0, 0.d0, 0.d0]
        ITA_lattice_vectors(2,1:3) = lattice_parameters(2) * [0.d0, 1.d0, 0.d0]                
        ITA_lattice_vectors(3,1:3) = lattice_parameters(3) * [0.d0, 0.d0, 1.d0]        
    end select

    if(debug) then
        write(6,'(a)') "toITAConventionalCell: ITA lattice vectors"
        write(6,'(3f16.9)') ITA_lattice_vectors(1,1:3)
        write(6,'(3f16.9)') ITA_lattice_vectors(2,1:3)
        write(6,'(3f16.9)') ITA_lattice_vectors(3,1:3)
    endif

end subroutine toITAConventionalCell

end module toITAConventionalCell_m