SUBROUTINE F_EQ
    USE VARIABLES
    
    
    !!$OMP PARALLEL DO PRIVATE(X,Y,Z,I,UXY,USQR) SHARED(WALL,UREAL,VREAL,WREAL,FEQ,WF,E,DXY)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(X,Y,Z,I,UXY,USQR) SHARED(WALL, UREAL, VREAL, WREAL, FEQ, WF, E, DXY) SCHEDULE(STATIC)
    DO Z=1,ZMAX
        DO Y=1,YMAX
            DO X=1,XMAX
                IF (WALL(X,Y,Z) .EQ. 0) THEN
                    
                    USQR(X,Y,Z) = U(X,Y,Z) * U(X,Y,Z) + V(X,Y,Z) * V(X,Y,Z) + W(X,Y,Z) * W(X,Y,Z)
    
                    
                    DO I=0,18
                        UXY(X,Y,Z) = U(X,Y,Z)*E(I,0) + V(X,Y,Z)*E(I,1) + W(X,Y,Z)*E(I,2)
                        FEQ(I,X,Y,Z) = WF(I) * DXY(X,Y,Z) * (1.0d0 + 3.0d0*UXY(X,Y,Z) + 4.5d0*UXY(X,Y,Z)*UXY(X,Y,Z) - 1.5d0*USQR(X,Y,Z))
                    END DO
                END IF
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE

    
