module symmetrizeLatticeVectors_m
    use constants, only: eps16
    use kinds, only: DP

    implicit none

    public :: symmetrizeLatticeVectors

    private

contains

subroutine symmetrizeLatticeVectors( &
    lattice_vectors, n_W, W, debug &
)
    integer, intent(in)     :: n_W, W(3,3,n_W)
    real(DP), intent(inout) :: lattice_vectors(3,3)
    logical, intent(in)     :: debug

    external :: dgesvd, dgemm
    integer  :: lwork, info
    real(DP) :: U(3,3), VT(3,3), S(3), R(3,3), P(3,3), temp_P(3,3), temp_S(3,3)
    real(DP), allocatable :: work(:)

    external :: dsyev
    integer  :: iw, ix, jx
    real(DP) :: G(3,3), temp_G(3,3), symm_G(3,3), theta(3), &
                Evec(3,3), Ev(3), temp_M(3,3), symm_P(3,3)

    real(DP) :: temp_lattice_vectors(3,3)

    do ix = 1, 3;do jx = 1, 3
        temp_lattice_vectors(ix,jx) = &
            lattice_vectors(jx,ix)
    enddo;enddo
    if(debug) then
        write(6,'(a)') "symmetrizeCrystal: Lattice vectors A"
        write(6,'(3f16.9)') temp_lattice_vectors(1,1:3)
        write(6,'(3f16.9)') temp_lattice_vectors(2,1:3)
        write(6,'(3f16.9)') temp_lattice_vectors(3,1:3)
    endif
    ! temp_lattice_vectors(1:3,1:3) = lattice_vectors(1:3,1:3)

    ! Polar decomposition of lattice vector matrix A = RP
    ! First use SVD to make A = USV^T
    ! Second R = UV^T, P = VSV^T
    lwork = 15
    allocate(work(lwork))
    call dgesvd('A','A',3,3,temp_lattice_vectors,3,S,U,3,VT,3,work,lwork,info)
    if(info .ne. 0) then
        write(6,'(a)') "symmetrizeLatticeVectors: SVD failed for lattice vector matrix"
        stop
    endif
    do ix = 1, 3;do jx = 1, 3
        R(ix,jx) = &
            U(ix,1) * VT(1,jx) + &
            U(ix,2) * VT(2,jx) + &
            U(ix,3) * VT(3,jx)
    enddo;enddo
    if(abs(R(1,1) - 1.d0) .gt. eps16) then
        theta(1) = acos(R(2,2))
    elseif(abs(R(2,2) - 1.d0) .gt. eps16) then
        theta(1) = acos(R(1,1))
    elseif(abs(R(3,3) - 1.d0) .gt. eps16) then
        theta(1) = acos(R(1,1))
    else
        write(6,'(a)') "symmetrizeLatticeVectors: R is not a rotation matrix"
        stop        
    endif
    if(debug) write(6,'(a,f16.9)') "symmetrizeLatticeVectors: theta (degree) = ", theta(1)*180.d0/(4.d0*atan(1.d0))
    temp_S(1,1:3) = [ S(1), 0.d0, 0.d0 ]
    temp_S(2,1:3) = [ 0.d0, S(2), 0.d0 ]
    temp_S(3,1:3) = [ 0.d0, 0.d0, S(3) ]
    write(6,*) temp_S
    write(6,*) 
    write(6,*) VT
    write(6,*) 
    do ix = 1, 3;do jx = 1, 3
        temp_P(ix,jx) = &
            VT(1,ix) * temp_S(1,jx) + &
            VT(2,ix) * temp_S(2,jx) + &
            VT(3,ix) * temp_S(3,jx)
    enddo;enddo
    write(6,*) temp_P
    do ix = 1, 3; do jx = 1, 3
        P(ix,jx) = &
            temp_P(ix,1) * VT(1,jx) + &
            temp_P(ix,2) * VT(2,jx) + &
            temp_P(ix,3) * VT(3,jx)
    enddo;enddo
    if(debug) then
        write(6,'(a)') "symmetrizeCrystal: R"
        write(6,'(3f16.9)') R(1,1:3)
        write(6,'(3f16.9)') R(2,1:3)
        write(6,'(3f16.9)') R(3,1:3)
        write(6,'(a)') "symmetrizeCrystal: P"
        write(6,'(3f16.9)') P(1,1:3)
        write(6,'(3f16.9)') P(2,1:3)
        write(6,'(3f16.9)') P(3,1:3)
        do ix = 1, 3;do jx = 1, 3
            temp_M(ix,jx) = &
                R(ix,1)*P(1,jx) + &
                R(ix,2)*P(2,jx) + &
                R(ix,3)*P(3,jx)
        enddo;enddo
        write(6,'(a)') "symmetrizeCrystal: RP"
        write(6,'(3f16.9)') temp_M(1,1:3)
        write(6,'(3f16.9)') temp_M(2,1:3)
        write(6,'(3f16.9)') temp_M(3,1:3)
    endif

    do ix = 1, 3;do jx = 1, 3
        temp_lattice_vectors(ix,jx) = &
            lattice_vectors(jx,ix)
    enddo;enddo

    do ix = 1, 3;do jx = 1, 3
        G(ix,jx) = &
            temp_lattice_vectors(1,ix) * temp_lattice_vectors(1,jx) + &
            temp_lattice_vectors(2,ix) * temp_lattice_vectors(2,jx) + &
            temp_lattice_vectors(3,ix) * temp_lattice_vectors(3,jx)
    enddo;enddo

    if(debug) then
        write(6,'(a)') "symmetrizeCrystal: Metric tensor G"
        write(6,'(3f16.9)') G(1,1:3)
        write(6,'(3f16.9)') G(2,1:3)
        write(6,'(3f16.9)') G(3,1:3)
    endif

    symm_G(1:3,1:3) = 0.d0
    do iw = 1, n_W
        
        do ix = 1, 3; do jx = 1, 3
            temp_G(ix,jx) = &
                W(1,ix,iw) * G(1,jx) + & 
                W(2,ix,iw) * G(2,jx) + & 
                W(3,ix,iw) * G(3,jx)
        enddo;enddo
        do ix = 1, 3; do jx = 1, 3
            symm_G(ix,jx) = &
                symm_G(ix,jx) + &
                temp_G(ix,1) * W(1,jx,iw) + &
                temp_G(ix,2) * W(2,jx,iw) + &
                temp_G(ix,3) * W(3,jx,iw)
        enddo;enddo
    enddo
    symm_G(1:3,1:3) = symm_G(1:3,1:3) / dble(n_W)
    if(debug) then
        write(6,'(a)') "symmetrizeLatticeVectors: Symmetrized Metric tensor "
        write(6,'(3f16.9)') symm_G(1,1:3)
        write(6,'(3f16.9)') symm_G(2,1:3)
        write(6,'(3f16.9)') symm_G(3,1:3)
    endif
    
    do ix = 1, 3;do jx = 1, 3
        lattice_vectors(ix,jx) = &
            R(ix,1) * symm_P(1,jx) + &
            R(ix,2) * symm_P(2,jx) + &
            R(ix,3) * symm_P(3,jx)
    enddo;enddo

    ! We want to finde the symmetrized P
    ! Diagonalize symmetrized G so that G = Udiag(x1,x2,x3)U^T
    ! A_symm = Udiag(sqrt(x1),sqrt(x2),sqrt(x3))
    Evec(1:3,1:3) = symm_G(1:3,1:3)
    lwork = -1
    deallocate(work)
    allocate(work(1))
    info = 1
    call dsyev('V','U',3,Evec,3,Ev,work,lwork,info)
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    call dsyev('V','U',3,Evec,3,Ev,work,lwork,info)
    if(info .ne. 0) then
        write(6,'(a)') "symmetrizeLatticeVectors: Diagonalization failed for symmetrized G"
        stop
    endif
    if(debug) then
        write(6,'(a)') "symmetrizeLatticeVectors: Diagonalization of symmetrized metric tensor"
        write(6,'(a)') "symmetrizeLatticeVectors: Eigenvectors"
        write(6,'(3f16.9)') Evec(1,1:3)
        write(6,'(3f16.9)') Evec(2,1:3)
        write(6,'(3f16.9)') Evec(3,1:3)
        write(6,'(a)') "symmetrizeLatticeVectors: Eigenvalues"
        write(6,'(3f16.9)') Ev(1:3)
    endif

    temp_M(1,1:3) = [ sqrt(Ev(1)),        0.d0,        0.d0 ]
    temp_M(2,1:3) = [        0.d0, sqrt(Ev(2)),        0.d0 ]
    temp_M(3,1:3) = [        0.d0,        0.d0, sqrt(Ev(3)) ]
    do ix = 1, 3;do jx = 1, 3
        temp_P(ix,jx) = &
            Evec(ix,1) * temp_M(1,jx) + &
            Evec(ix,2) * temp_M(2,jx) + &
            Evec(ix,3) * temp_M(3,jx)
    enddo;enddo
    do ix = 1, 3;do jx = 1, 3
        symm_P(ix,jx) = &
            temp_P(ix,1) * Evec(jx,1) + &
            temp_P(ix,2) * Evec(jx,2) + &
            temp_P(ix,3) * Evec(jx,3)
    enddo;enddo
    if(debug) then
        write(6,'(a)') "symmetrizeLatticeVectors: Symmetrized P"
        write(6,'(3f16.9)') symm_P(1,1:3)
        write(6,'(3f16.9)') symm_P(2,1:3)
        write(6,'(3f16.9)') symm_P(3,1:3)
    endif

    if(debug) then
        write(6,'(a)') "symmetrizeLatticeVectors: Symmetrized Lattice vector "
        write(6,'(3f16.9)') lattice_vectors(1:3,1)
        write(6,'(3f16.9)') lattice_vectors(1:3,2)
        write(6,'(3f16.9)') lattice_vectors(1:3,3)
    endif

    do ix = 1, 3;do jx = 1, 3
        temp_lattice_vectors(ix,jx) = &
            lattice_vectors(jx,ix)
    enddo;enddo
    lattice_vectors = temp_lattice_vectors

end subroutine symmetrizeLatticeVectors

end module symmetrizeLatticeVectors_m