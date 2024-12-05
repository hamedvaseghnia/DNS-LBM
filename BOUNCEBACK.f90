SUBROUTINE BOUNCEBACK
    USE VARIABLES
    IMPLICIT NONE
    
    ! Parallelize the outermost loop (Z loop)
    !!$OMP PARALLEL DO PRIVATE(X,Y,Z,I) SHARED(FIN,FOUT,WALL,OPPOSITE)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(X,Y,Z,I) SHARED(FIN,FOUT,WALL,OPPOSITE) SCHEDULE(STATIC)
    DO Z=1,ZMAX
        DO Y=1,YMAX
            DO X=1,XMAX
                IF (WALL(X,Y,Z).EQ. 1) THEN   
                    ! Loop over the distribution function directions
                    DO I=0,18
                        FOUT(I,X,Y,Z) = FIN(OPPOSITE(I),X,Y,Z)
                    END DO
                END IF
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE
    