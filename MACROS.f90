SUBROUTINE MACROS
    USE VARIABLES
    
    !!$OMP PARALLEL DO PRIVATE(X,Y,Z) SHARED(FIN,DXY,U,V,W,WALL)
     !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(X,Y,Z) SHARED(FIN,DXY,U,V,W,WALL) SCHEDULE(STATIC)
    DO Z=1,ZMAX
        DO Y = 1,YMAX
            DO X = 1,XMAX
                
                IF (WALL(X,Y,Z) .EQ. 0) THEN
                    DXY(X,Y,Z) = FIN(0,X,Y,Z) + FIN(1,X,Y,Z) + FIN(2,X,Y,Z) + FIN(3,X,Y,Z) + FIN(4,X,Y,Z) + FIN(5,X,Y,Z) + FIN(6,X,Y,Z) + FIN(7,X,Y,Z) + FIN(8,X,Y,Z) + FIN(9,X,Y,Z) + FIN(10,X,Y,Z) + FIN(11,X,Y,Z) + FIN(12,X,Y,Z) + FIN(13,X,Y,Z) + FIN(14,X,Y,Z) + FIN(15,X,Y,Z) + FIN(16,X,Y,Z) + FIN(17,X,Y,Z) + FIN(18,X,Y,Z)
    
                    U(X,Y,Z) = (FIN(1,X,Y,Z) - FIN(2,X,Y,Z) + FIN(7,X,Y,Z) - FIN(8,X,Y,Z) + FIN(9,X,Y,Z) - FIN(10,X,Y,Z) + FIN(13,X,Y,Z) - FIN(14,X,Y,Z) + FIN(15,X,Y,Z) - FIN(16,X,Y,Z)) / DXY(X,Y,Z)
    
                    V(X,Y,Z) = (FIN(3,X,Y,Z) - FIN(4,X,Y,Z) + FIN(7,X,Y,Z) - FIN(8,X,Y,Z) + FIN(11,X,Y,Z) - FIN(12,X,Y,Z) - FIN(13,X,Y,Z) + FIN(14,X,Y,Z) + FIN(17,X,Y,Z) - FIN(18,X,Y,Z)) / DXY(X,Y,Z)
    
                    W(X,Y,Z) = (FIN(5,X,Y,Z) - FIN(6,X,Y,Z) + FIN(9,X,Y,Z) - FIN(10,X,Y,Z) + FIN(11,X,Y,Z) - FIN(12,X,Y,Z) - FIN(15,X,Y,Z) + FIN(16,X,Y,Z) - FIN(17,X,Y,Z) + FIN(18,X,Y,Z)) / DXY(X,Y,Z)
                END IF
            END DO
        END DO
    END DO
    !!$OMP END PARALLEL DO
    !$OMP END PARALLEL DO 
    
    END SUBROUTINE
    
