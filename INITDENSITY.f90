SUBROUTINE INITDENSITY
USE VARIABLES
IMPLICIT NONE

DO Z=1,ZMAX
    DO Y=1,YMAX
         DO X=1,XMAX
         USQR(X,Y,Z) = U(X,Y,Z) * U(X,Y,Z) + V(X,Y,Z) * V(X,Y,Z) + W(X,Y,Z)*W(X,Y,Z)
                DO I=0,18
                UXY(X,Y,Z)=  U(X,Y,Z)*E(I,0)  +  V(X,Y,Z)*E(I,1)  +  W(X,Y,Z)*E(I,2)
                FIN(I,X,Y,Z)= WF(I)*DXY(X,Y,Z)*(1.0D0+3.0D0*UXY(X,Y,Z)+4.5D0*UXY(X,Y,Z)*UXY(X,Y,Z)-1.5D0*USQR(X,Y,Z))
                END DO     
        END DO
    END DO
END DO

print*,'Initialization:  DONE'
END SUBROUTINE