SUBROUTINE WALL_COORDINATES
  USE VARIABLES
  IMPLICIT NONE
  
  DO Z=1,ZMAX
      DO Y=1,YMAX
          DO X=1,XMAX
            WALL(X,Y,Z)=0
          END DO
      END DO
  END DO
  
  
  
     print*,''
      PRINT*,"Total number of nodes=",xmax*ymax*zmax
  
  END SUBROUTINE