SUBROUTINE PRESS_BOUND
    USE VARIABLES
    IMPLICIT NONE
    
    !=========================XMIN BOUNDARY CONDITION============================
    
    X = 1
    ! Parallelize the outermost loop over Z for the XMIN boundary
    !!$OMP PARALLEL DO PRIVATE(Y,Z) SHARED(X, WALL, FIN, FOUT, UXY, DXY)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(Y,Z) SHARED(X, WALL, FIN, FOUT, UXY, DXY) SCHEDULE(STATIC)
    DO Z = 1, ZMAX
        DO Y = 1, YMAX
            IF (WALL(X, Y, Z) .EQ. 0) THEN
               !IF (STEP .LE.10000)THEN
                !UXY(X,Y,Z) = 1 * 0.0311111111111111d0 * (1.0d0 - 4.0d0 * ( ((Z - (75.5d0 / 2.0d0 + 1.0d0))**2 + (Y - (75.5d0 / 2.0d0 + 1.0d0))**2) ) / 75.5d0**2) 
               !ELSE IF (STEP .LE. 20000) THEN
                   !UXY(X,Y,Z) = 2 * 0.0311111111111111d0*STEP/20000 * (1.0d0 - 4.0d0 * ( ((Z - (75.5d0 / 2.0d0 + 1.0d0))**2 + (Y - (75.5d0 / 2.0d0 + 1.0d0))**2) ) / 75.5d0**2)     
               !ELSE
                UXY(X,Y,Z) = 2 * 0.0288888888888889d0 * (1.0d0 - 4.0d0 * ( ((Z - (75.5d0 / 2.0d0 + 1.5d0))**2 + (Y - (75.5d0 / 2.0d0 + 1.5d0))**2) ) / 75.5d0**2)
               !END IF     
                   DXY(X, Y, Z) = (FIN(0, X, Y, Z) + FIN(3, X, Y, Z) + FIN(4, X, Y, Z) + FIN(5, X, Y, Z) + FIN(6, X, Y, Z) + FIN(18, X, Y, Z) + FIN(11, X, Y, Z) + FIN(17, X, Y, Z) + FIN(12, X, Y, Z) + 2.D0 * (FIN(2, X, Y, Z) + FIN(14, X, Y, Z) + FIN(8, X, Y, Z) + FIN(16, X, Y, Z) + FIN(10, X, Y, Z))) / (1D0 + UXY(X, Y, Z))   
    
                FIN(1, X, Y, Z) = FIN(2, X, Y, Z) + (1.D0 / 3.D0) * UXY(X, Y, Z) * DXY(X, Y, Z)
                FIN(13, X, Y, Z) = FIN(14, X, Y, Z) + (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) + (1.D0 / 2.D0) * (FIN(3, X, Y, Z) + FIN(11, X, Y, Z) + FIN(17, X, Y, Z) - FIN(4, X, Y, Z) - FIN(18, X, Y, Z) - FIN(12, X, Y, Z))
                FIN(7, X, Y, Z) = FIN(8, X, Y, Z) + (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) - (1.D0 / 2.D0) * (FIN(3, X, Y, Z) + FIN(11, X, Y, Z) + FIN(17, X, Y, Z) - FIN(4, X, Y, Z) - FIN(18, X, Y, Z) - FIN(12, X, Y, Z))
                FIN(9, X, Y, Z) = FIN(10, X, Y, Z) + (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) - (1.D0 / 2.D0) * (FIN(5, X, Y, Z) + FIN(18, X, Y, Z) + FIN(11, X, Y, Z) - FIN(6, X, Y, Z) - FIN(17, X, Y, Z) - FIN(12, X, Y, Z))
                FIN(15, X, Y, Z) = FIN(16, X, Y, Z) + (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) + (1.D0 / 2.D0) * (FIN(5, X, Y, Z) + FIN(18, X, Y, Z) + FIN(11, X, Y, Z) - FIN(6, X, Y, Z) - FIN(17, X, Y, Z) - FIN(12, X, Y, Z))
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    !=========================XMAX BOUNDARY CONDITION============================
    
    X = XMAX
    ! Parallelize the outermost loop over Z for the XMAX boundary
    !!$OMP PARALLEL DO PRIVATE(Y,Z) SHARED(X, WALL, FIN, FOUT, UXY, DXY)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(Y,Z) SHARED(X, WALL, FIN, FOUT, UXY, DXY) SCHEDULE(STATIC)
    DO Z = 1, ZMAX
        DO Y = 1, YMAX
            IF (WALL(X, Y, Z) .EQ. 0) THEN
                DXY(X, Y, Z) = 1D0
                UXY(X, Y, Z) = -1.D0 + (FIN(0, X, Y, Z) + FIN(3, X, Y, Z) + FIN(4, X, Y, Z) + FIN(5, X, Y, Z) + FIN(6, X, Y, Z) + FIN(18, X, Y, Z) + FIN(11, X, Y, Z) + FIN(17, X, Y, Z) + FIN(12, X, Y, Z) + 2.D0 * (FIN(2, X, Y, Z) + FIN(14, X, Y, Z) + FIN(8, X, Y, Z) + FIN(16, X, Y, Z) + FIN(10, X, Y, Z))) / DXY(X, Y, Z)
    
                FIN(2, X, Y, Z)  = FIN(1, X, Y, Z) - (1.D0 / 3.D0) * UXY(X, Y, Z) * DXY(X, Y, Z)
                FIN(14, X, Y, Z) = FIN(13, X, Y, Z) - (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) - (1.D0 / 2.D0) * (FIN(3, X, Y, Z) + FIN(11, X, Y, Z) + FIN(17, X, Y, Z) - FIN(4, X, Y, Z) - FIN(18, X, Y, Z) - FIN(12, X, Y, Z))    
                FIN(8, X, Y, Z)  = FIN(7, X, Y, Z) - (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) + (1.D0 / 2.D0) * (FIN(3, X, Y, Z) + FIN(11, X, Y, Z) + FIN(17, X, Y, Z) - FIN(4, X, Y, Z) - FIN(18, X, Y, Z) - FIN(12, X, Y, Z))    
                FIN(10, X, Y, Z) = FIN(9, X, Y, Z) - (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) + (1.D0 / 2.D0) * (FIN(5, X, Y, Z) + FIN(18, X, Y, Z) + FIN(11, X, Y, Z) - FIN(6, X, Y, Z) - FIN(17, X, Y, Z) - FIN(12, X, Y, Z))    
                FIN(16, X, Y, Z) = FIN(15, X, Y, Z) - (1.D0 / 6.D0) * UXY(X, Y, Z) * DXY(X, Y, Z) - (1.D0 / 2.D0) * (FIN(5, X, Y, Z) + FIN(18, X, Y, Z) + FIN(11, X, Y, Z) - FIN(6, X, Y, Z) - FIN(17, X, Y, Z) - FIN(12, X, Y, Z))     
            END IF
        END DO
    END DO

    !FIN(2, X-1, Y, Z)!
    !FIN(14,X-1, Y, Z)!
    !FIN(8, X-1, Y, Z)!
    !FIN(10,X-1, Y, Z)!
    !FIN(16,X-1, Y, Z)!


    !$OMP END PARALLEL DO
    
    END SUBROUTINE
    