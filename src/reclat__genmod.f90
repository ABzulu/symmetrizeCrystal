        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  3 18:51:32 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECLAT__genmod
          INTERFACE 
            SUBROUTINE RECLAT(A,B,IOPT)
              REAL(KIND=8), INTENT(IN) :: A(3,3)
              REAL(KIND=8), INTENT(OUT) :: B(3,3)
              INTEGER(KIND=4), INTENT(OUT) :: IOPT
            END SUBROUTINE RECLAT
          END INTERFACE 
        END MODULE RECLAT__genmod
