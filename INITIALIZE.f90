SUBROUTINE INITIALIZE
  USE VARIABLES
  IMPLICIT NONE
    
      
  DO Z=1,ZMAX
    DO Y=1,YMAX
    DO X=1,XMAX 
          T=TR*TC
          IF (WALL (X,Y,Z) .EQ. 0 ) THEN 
  
          
            U(X,Y,Z)= (0.025d0)*dsin((2d0*pi/ymax)*x)*dcos((2d0*pi/ymax)*y)*dcos((2d0*pi/ymax)*z)
            V(X,Y,Z)=-(0.025d0)*dcos((2d0*pi/ymax)*x)*dsin((2d0*pi/ymax)*y)*cos((2d0*pi/ymax)*z)
            W(X,Y,Z)= 0d0
          
          dxy(x,y,z) = 1d0!((1d0/3d0) + ((0.1d0*CS)**2d0/16d0) * &
               !(dcos((4d0 * pi / ymax) * x) + dcos((4d0 * pi / ymax) * y)) * &
               !(dcos((4d0 * pi / ymax) * z) + 2d0))/(CS**2d0)
          
          RELAXATION(X,Y,Z)= 3d0*0.00127323954473516d0 + 0.5d0  !0.520750565d0
  
          END IF
      END DO
    END DO
  END DO
  
  
!  DO Z=1,ZMAX
!    DO Y=1,YMAX
!        DO X=1,XMAX
!        
!
!            YT = MOD(Y, YMAX) + 1 
!            XR = MOD(X, XMAX) + 1 
!            YB = YMAX - MOD(YMAX + 1 - Y, YMAX) 
!            XL = XMAX - MOD(XMAX + 1 - X, XMAX)
!            ZI = MOD(Z, ZMAX) + 1 
!            ZO = ZMAX - MOD(ZMAX + 1 - Z, ZMAX)
!
!
!            s_11(X,Y,Z) = (U(XR,Y,Z) - U(XL,Y,Z)) / 2D0
!            s_12(X,Y,Z) = (((U(X,YT,Z) - U(X,YB,Z)) / 2D0) + ((V(XL,Y,Z) - V(XR,Y,Z)) / 2D0)) / 2D0
!            s_13(X,Y,Z) = (((U(X,Y,ZI) - U(X,YB,ZO)) / 2D0) + ((W(XL,Y,Z) - W(XR,Y,Z)) / 2D0)) / 2D0
!            s_22(X,Y,Z) = (V(X,YT,Z) - V(X,YB,Z)) / 2D0
!            s_23(X,Y,Z) = (((V(X,Y,ZI) - V(X,Y,ZO)) / 2D0) + ((W(X,YT,Z) - W(X,YB,Z)) / 2D0)) / 2D0
!            s_33(X,Y,Z) = (W(X,Y,ZI) - W(X,Y,ZO)) / 2D0
!
!        END DO
!    END DO
!END DO
!
!
!DO Z=1,ZMAX
!    DO Y=1,YMAX
!        DO X=1,XMAX
!            S_squared(X,Y,Z) = s_11(X,Y,Z)**2D0 + s_22(X,Y,Z)**2D0 + s_33(X,Y,Z)**2D0 + 2D0 * (s_12(X,Y,Z)**2D0 + s_13(X,Y,Z)**2D0 + s_23(X,Y,Z)**2D0)
!        END DO
!    END DO
!END DO
!
!
!DO Z=1,ZMAX
!    DO Y=1,YMAX
!        DO X=1,XMAX
!            omega_squared(X,Y,Z) = 2D0 * (((U(X,YT,Z) - U(X,YB,Z)) / 2D0 - (V(XL,Y,Z) - V(XR,Y,Z)) / 2D0)**2D0  + &
!                                         ((U(X,Y,ZI) - U(X,Y,ZO)) / 2D0 - (W(XL,Y,Z) - W(XR,Y,Z)) / 2D0)**2D0  + &
!                                         ((V(X,Y,ZI) - V(X,Y,ZO)) / 2D0 - (W(X,YT,Z) - W(X,YB,Z)) / 2D0)**2D0)
!        END DO
!    END DO
!END DO
!
!
!DO Z=1,ZMAX
!    DO Y=1,YMAX
!        DO X=1,XMAX
!            Q_crit(X,Y,Z) = 0.5D0 * (omega_squared(X,Y,Z) - S_squared(X,Y,Z))
!        END DO
!    END DO
!END DO
!
!  
!  
!  OPEN(20,FILE="INTITAL.PLT")
!  
!  WRITE(20,*)"VARIABLES= X,Y,Z,DXY,Q"
!  WRITE(20,*)"ZONE I=" ,XMAX , "J=",YMAX , "K=",ZMAX , "F=POINT"
!  
!  DO Z=1,ZMAX
!      DO Y=1,YMAX
!           DO X=1,XMAX
!           WRITE(20,*)X,Y,Z,DXY(X,Y,Z),Q_crit(X,Y,Z)
!           END DO
!      END DO
!  END DO         
!  
!  
!  
!  DO Z=1,ZMAX
!      DO Y=1,YMAX
!           DO X=1,XMAX
!                  SS(1,X,Y,Z) =   0.D0    !rho
!                  SS(2,X,Y,Z) =   1.19D0    !e
!                  SS(3,X,Y,Z) =   1.2D0    !epsilon 
!                  SS(4,X,Y,Z) =   0.D0    !jx
!                  SS(5,X,Y,Z) =   1.2D0    !qx
!                  SS(6,X,Y,Z) =   0.D0    !jy
!                  SS(7,X,Y,Z) =   1.2D0    !qy
!                  SS(8,X,Y,Z) =   0.D0    !jz 
!                  SS(9,X,Y,Z) =   1.2D0    !qz 
!                  SS(10,X,Y,Z)=  1.D0/RELAXATION(X,Y,Z)    !pxx
!                  SS(11,X,Y,Z)=  1.4D0    !PIxx 
!                  SS(12,X,Y,Z)=  1.D0/RELAXATION(X,Y,Z)    !pww
!                  SS(13,X,Y,Z)=  1.4D0    !PIww
!                  SS(14,X,Y,Z)=  1.D0/RELAXATION(X,Y,Z)    !pxz
!                  SS(15,X,Y,Z)=  1.D0/RELAXATION(X,Y,Z)    !pyz
!                  SS(16,X,Y,Z)=  1.D0/RELAXATION(X,Y,Z) 
!                  SS(17,X,Y,Z)=  1.6D0   !mx
!                  SS(18,X,Y,Z)=  1.6D0   !my
!                  SS(19,X,Y,Z)=  1.6D0   !mz
!           END DO
!      END DO
!  END DO  
  
  
  
  
  
  
  !OPEN(3,FILE="DELTA P.PLT")
  !WRITE(3,*)"VARIABLES= ITERATION,DELTA_P"
  
  END SUBROUTINE