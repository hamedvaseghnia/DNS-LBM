SUBROUTINE COLLISION
    USE VARIABLES
    
    !! Initialize arrays to zero
    !MFIN(:,:,:,:)=0.D0
    !MFEQ(:,:,:,:)=0.D0
    !TRANSFORMED(:,:,:,:)=0.D0
    !
    !Sab_XX(:,:,:)=0D0      
    !Sab_YY(:,:,:)=0D0      
    !Sab_XY(:,:,:)=0D0
    !Sab_XZ(:,:,:)=0D0
    !Sab_YZ(:,:,:)=0D0
    !Sab_ZZ(:,:,:)=0D0
    !
    !!$OMP PARALLEL DO PRIVATE(Z, Y, X, I, II, JJ, TMP_RST, NEQ, SHEAR, TMP_RELAX, N_INDEX, INVARIANT) SHARED(WALL, FIN, FEQ, E, DXY, RELAXATION, M, MINV, SS, FOUT, MRT_FX, MRT_FY, MRT_FZ, MFIN, MFEQ, TRANSFORMED, Sab_XX, Sab_YY, Sab_XY, Sab_XZ, Sab_YZ, Sab_ZZ, RST_XX, RST_YY, RST_XY, RST_XZ, RST_YZ, RST_ZZ, VISCOSITY)
    !DO Z=1,ZMAX
    !    DO Y=1,YMAX
    !        DO X=1,XMAX
    !            IF (WALL(X,Y,Z) .EQ. 0) THEN  
    !                DO I=0,18
    !                    NEQ(I,X,Y,Z) = FIN(I,X,Y,Z) - FEQ(I,X,Y,Z)
    !                    Sab_XX(X,Y,Z) = Sab_XX(X,Y,Z) + NEQ(I,X,Y,Z) * E(I,0) * E(I,0)       
    !                    Sab_YY(X,Y,Z) = Sab_YY(X,Y,Z) + NEQ(I,X,Y,Z) * E(I,1) * E(I,1)  
    !                    Sab_ZZ(X,Y,Z) = Sab_ZZ(X,Y,Z) + NEQ(I,X,Y,Z) * E(I,2) * E(I,2)
    !                    Sab_XY(X,Y,Z) = Sab_XY(X,Y,Z) + NEQ(I,X,Y,Z) * E(I,0) * E(I,1)
    !                    Sab_XZ(X,Y,Z) = Sab_XZ(X,Y,Z) + NEQ(I,X,Y,Z) * E(I,0) * E(I,2)
    !                    Sab_YZ(X,Y,Z) = Sab_YZ(X,Y,Z) + NEQ(I,X,Y,Z) * E(I,1) * E(I,2)
    !                END DO
    !
    !                TMP_RST(X,Y,Z) = -1.5D0 / (DXY(X,Y,Z) * RELAXATION(X,Y,Z))
    !                RST_XX(X,Y,Z) = TMP_RST(X,Y,Z) * Sab_XX(X,Y,Z)
    !                RST_YY(X,Y,Z) = TMP_RST(X,Y,Z) * Sab_YY(X,Y,Z)
    !                RST_XY(X,Y,Z) = TMP_RST(X,Y,Z) * Sab_XY(X,Y,Z)
    !                RST_XZ(X,Y,Z) = TMP_RST(X,Y,Z) * Sab_XZ(X,Y,Z)
    !                RST_YZ(X,Y,Z) = TMP_RST(X,Y,Z) * Sab_YZ(X,Y,Z)
    !                RST_ZZ(X,Y,Z) = TMP_RST(X,Y,Z) * Sab_ZZ(X,Y,Z)   
    !
    !                INVARIANT(X,Y,Z) = RST_XX(X,Y,Z)**2 + 2D0*RST_XY(X,Y,Z)**2 + RST_YY(X,Y,Z)**2 + 2D0*RST_XZ(X,Y,Z)**2 + 2D0*RST_YZ(X,Y,Z)**2 + RST_ZZ(X,Y,Z)**2
    !                SHEAR(X,Y,Z) = DSQRT(2D0 * INVARIANT(X,Y,Z))
    !                N_INDEX = 0.8D0
    !
    !                ! Ensure shear is not too small
    !                !IF (SHEAR(X,Y,Z) .LT. 1D-6) THEN
    !                !    SHEAR(X,Y,Z) = 1D-6
    !                !END IF
    !                VISCOSITY(X,Y,Z) = 0.00127323954473516d0!0.1018591*YMAX/(600d0*pi)  !0.00691685515259379d0!0.01D0 !* ABS(SHEAR(X,Y,Z)**(N_INDEX - 1D0))
    !                RELAXATION(X,Y,Z) = (6D0 * VISCOSITY(X,Y,Z) + 1D0) / 2D0  
    !            END IF
    !
    !            TMP_RELAX = 1.D0 / RELAXATION(X,Y,Z) 
    !            SS(10,X,Y,Z) = TMP_RELAX    
    !            SS(12,X,Y,Z) = TMP_RELAX    
    !            SS(14,X,Y,Z) = TMP_RELAX    
    !            SS(15,X,Y,Z) = TMP_RELAX    
    !            SS(16,X,Y,Z) = TMP_RELAX
    !
    !            ! Matrix operations
    !            DO II=1,19
    !                DO JJ=1,19
    !                    MFIN(II-1,X,Y,Z) = MFIN(II-1,X,Y,Z) + M(II,JJ) * FIN(JJ-1,X,Y,Z)
    !                    MFEQ(II-1,X,Y,Z) = MFEQ(II-1,X,Y,Z) + M(II,JJ) * FEQ(JJ-1,X,Y,Z)
    !                END DO
    !            END DO
    !
    !            DO II=1,19
    !                RELAXEDM(II-1,X,Y,Z) = SS(II,X,Y,Z) * (MFIN(II-1,X,Y,Z) - MFEQ(II-1,X,Y,Z))
    !            END DO
    !
    !            DO II=1,19
    !                DO JJ=1,19
    !                    TRANSFORMED(II-1,X,Y,Z) = TRANSFORMED(II-1,X,Y,Z) + MINV(II,JJ) * RELAXEDM(JJ-1,X,Y,Z)
    !                END DO
    !            END DO
    !
    !            IF (WALL(X,Y,Z) .EQ. 0) THEN 
    !                DO I=0,18
    !                    FOUT(I,X,Y,Z) = FIN(I,X,Y,Z) - TRANSFORMED(I,X,Y,Z) + (MRT_FX(I,X,Y,Z) + MRT_FY(I,X,Y,Z) + MRT_FZ(I,X,Y,Z))
    !                END DO
    !            END IF
    !        END DO
    !    END DO
    !END DO
    !!$OMP END PARALLEL DO
    !
    !END SUBROUTINE
    
     !========================== SRT IMPLEMENTATIION=====================================
    
    
    !$OMP PARALLEL DO PRIVATE(Z, Y, X, I) SHARED(WALL, FIN, FEQ, FOUT, RELAXATION)
    DO Z = 1, ZMAX
        DO Y = 1, YMAX
            DO X = 1, XMAX
                IF (WALL(X,Y,Z) .EQ. 0) THEN
                    VISCOSITY(X,Y,Z) = 0.00127323954473516d0!0.1018591*YMAX/(600d0*pi)  !0.00691685515259379d0!0.01D0 !* ABS(SHEAR(X,Y,Z)**(N_INDEX - 1D0))
                    RELAXATION(X,Y,Z) = (6D0 * VISCOSITY(X,Y,Z) + 1D0) / 2D0 
                    DO I = 0, 18
                        FOUT(I,X,Y,Z) = FIN(I,X,Y,Z) - (1d0 / RELAXATION(X,Y,Z)) * (FIN(I,X,Y,Z) - FEQ(I,X,Y,Z))
                    END DO
                END IF
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE
    