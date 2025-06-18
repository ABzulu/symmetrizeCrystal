module findRotationalSymmetry_m
    use constants, only: eps16

    implicit none

    public :: findRotationalSymmetry

    private

contains

subroutine findRotationalSymmetry( &
    lattice_vectors, lattice_tol, W, n_W, debug &
)
    double precision, intent(in) :: lattice_vectors(3,3), lattice_tol
    integer, intent(out) :: W(3,3,48)
    integer, intent(out) :: n_W
    logical, intent(in) :: debug

    integer :: ROT(3,3), det

    integer :: icounter, ix, jx
    double precision :: &
        G(3,3), temp_G(3,3), G_tilda(3,3), diag(3), off_diag(6), &
        cos1, cos2, theta
    integer :: i11, i12, i13, i21, i22, i23, i31, i32, i33, i, j

    do ix = 1, 3;do jx = 1, 3
        G(ix,jx) = &
            lattice_vectors(ix,1) * lattice_vectors(jx,1) + &
            lattice_vectors(ix,2) * lattice_vectors(jx,2) + &
            lattice_vectors(ix,3) * lattice_vectors(jx,3)
    enddo;enddo

    W = 0
    n_W = 0
    ! if(debug) write(6,'(a)') "findRotationalSymmetry: Candidate W"
    do i11 = -2,2; do i12 = -2,2; do i13 = -2,2
    do i21 = -2,2; do i22 = -2,2; do i23 = -2,2
    do i31 = -2,2; do i32 = -2,2; do i33 = -2,2
        det = &
            i11*(i22*i33 - i23*i32) - &
            i12*(i21*i33 - i23*i31) + &
            i13*(i21*i32 - i22*i31)
        if((.not.(det .eq. 1)) .and. (.not.(det .eq. -1))) cycle

        ROT(1,1:3) = [i11, i12, i13]
        ROT(2,1:3) = [i21, i22, i23]
        ROT(3,1:3) = [i31, i32, i33]

        ! if(debug) then
        !     write(6,'(3i4)') ROT(1,1:3)
        !     write(6,'(3i4)') ROT(2,1:3)
        !     write(6,'(3i4)') ROT(3,1:3)
        ! endif

        do i = 1, 3
            do j = 1, 3
                temp_G(i,j) = &
                    ROT(1,i) * G(1,j) + &
                    ROT(2,i) * G(2,j) + &
                    ROT(3,i) * G(3,j)
            enddo
        enddo
        do i = 1, 3
            do j = 1, 3
                G_tilda(i,j) = &
                    temp_G(i,1) * ROT(1,j) + &
                    temp_G(i,2) * ROT(2,j) + &
                    temp_G(i,3) * ROT(3,j)
            enddo
        enddo

        ! Diagonal elements represent lengths
        do i = 1, 3
            diag(i) = abs(G_tilda(i,i) - G(i,i))
        enddo
        ! Off-diagonal elements represent angles
        ! https://arxiv.org/html/1808.01590v2 Section V.1 Step(f)
        icounter = 0
        do i = 1, 3
            do j = 1, 3
                if(i == j) cycle
                icounter = icounter + 1
                cos1 = G_tilda(i,j)/sqrt(G_tilda(i,i)*G_tilda(j,j))
                cos2 = G(i,j)/sqrt(G(i,i)*G(j,j))
                theta = &
                    acos(max(-1.d0,min(1.d0,cos1))) - &
                    acos(max(-1.d0,min(1.d0,cos2)))
                off_diag(icounter) = &
                    sin(abs(theta))*sqrt(0.25d0*(G_tilda(i,i)+G(i,i))*(G_tilda(j,j)+G(j,j)))
            enddo
        enddo
        ! if(debug) write(6,*) "findRotationalSymmetry: diag, off_diag"
        ! if(debug) write(6,'(3f16.9)') diag(1:3)
        ! if(debug) write(6,'(3f16.9)') off_diag(1:3)
        ! if(debug) write(6,'(3f16.9)') off_diag(4:6)
        if(((maxval(diag)/sqrt(G(3,3))) .gt. lattice_tol) .or. &
           ((maxval(off_diag)/sqrt(G(3,3))) .gt. lattice_tol)) cycle

        n_W = n_W + 1
        if(n_W .gt. 48) then
            write(6,*) "findRotationalSymmetry: Too many symmetry operations"
            stop
        endif
        W(1:3,1:3,n_W) = ROT(1:3,1:3)

        if(debug) then
            write(6,'(a)') "findRotationalSymmetry: Found symmetric rotational operators"
            write(6,'(a,i4)') "findRotationalSymmetry: Rotational operator index = ", n_W
            write(6,'(3i6)') W(1,1:3,n_W)
            write(6,'(3i6)') W(2,1:3,n_W)
            write(6,'(3i6)') W(3,1:3,n_W)
        endif
    enddo;enddo;enddo;enddo;enddo;enddo;enddo;enddo;enddo

end subroutine findRotationalSymmetry

end module findRotationalSymmetry_m