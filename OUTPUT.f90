SUBROUTINE OUTPUT
    USE VARIABLES
    IMPLICIT NONE
    
    
    
    
    
    !DO Z=1,ZMAX
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
    !            s_11(X,Y,Z) = (UREAL(XR,Y,Z) - UREAL(XL,Y,Z)) / 2D0
    !            s_12(X,Y,Z) = (((UREAL(X,YT,Z) - UREAL(X,YB,Z)) / 2D0) + ((VREAL(XL,Y,Z) - VREAL(XR,Y,Z)) / 2D0)) / 2D0
    !            s_13(X,Y,Z) = (((UREAL(X,Y,ZI) - UREAL(X,YB,ZO)) / 2D0) + ((WREAL(XL,Y,Z) - WREAL(XR,Y,Z)) / 2D0)) / 2D0
    !            s_22(X,Y,Z) = (VREAL(X,YT,Z) - VREAL(X,YB,Z)) / 2D0
    !            s_23(X,Y,Z) = (((VREAL(X,Y,ZI) - VREAL(X,Y,ZO)) / 2D0) + ((WREAL(X,YT,Z) - WREAL(X,YB,Z)) / 2D0)) / 2D0
    !            s_33(X,Y,Z) = (WREAL(X,Y,ZI) - WREAL(X,Y,ZO)) / 2D0
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
    !            omega_squared(X,Y,Z) = 2D0 * (((UREAL(X,YT,Z) - UREAL(X,YB,Z)) / 2D0 - (VREAL(XL,Y,Z) - VREAL(XR,Y,Z)) / 2D0)**2D0  + &
    !                                         ((UREAL(X,Y,ZI) - UREAL(X,Y,ZO)) / 2D0 - (WREAL(XL,Y,Z) - WREAL(XR,Y,Z)) / 2D0)**2D0  + &
    !                                         ((VREAL(X,Y,ZI) - VREAL(X,Y,ZO)) / 2D0 - (WREAL(X,YT,Z) - WREAL(X,YB,Z)) / 2D0)**2D0)
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
    
    
    
                ! Periodic boundary conditions
                !YT = MOD(Y, YMAX) + 1
                !XR = MOD(X, XMAX) + 1
                !YB = MOD(YMAX + Y - 2, YMAX) + 1
                !XL = MOD(XMAX + X - 2, XMAX) + 1
                !ZI = MOD(Z, ZMAX) + 1
                !ZO = MOD(ZMAX + Z - 2, ZMAX) + 1
                
                ! Kinetic energy dissipation rate calculation
                !local_energy(X,Y,Z) = 2D0 * 0.01D0 * ( &
                !    ((UREAL(XR,Y,Z) - UREAL(XL,Y,Z)) / 2D0)**2D0 + &  ! dUx/dx
                !    ((UREAL(X,YT,Z) - UREAL(X,YB,Z)) / 2D0)**2D0 + &  ! dUx/dy
                !    ((UREAL(X,Y,ZI) - UREAL(X,Y,ZO)) / 2D0)**2D0 + &  ! dUx/dz
                !    ((VREAL(XR,Y,Z) - VREAL(XL,Y,Z)) / 2D0)**2D0 + &  ! dVy/dx
                !    ((VREAL(X,YT,Z) - VREAL(X,YB,Z)) / 2D0)**2D0 + &  ! dVy/dy
                !    ((VREAL(X,Y,ZI) - VREAL(X,Y,ZO)) / 2D0)**2D0 + &  ! dVy/dz
                !    ((WREAL(XR,Y,Z) - WREAL(XL,Y,Z)) / 2D0)**2D0 + &  ! dWz/dx
                !    ((WREAL(X,YT,Z) - WREAL(X,YB,Z)) / 2D0)**2D0 + &  ! dWz/dy
                !    ((WREAL(X,Y,ZI) - WREAL(X,Y,ZO)) / 2D0)**2D0 )    ! dWz/dz
    
    
    
    
    
    
    
    
    
    
    
    PMAX=MAXVAL(U(:,:,:))
    DP=PMAX
    
    UMAX=DSQRT(MAXVAL(USQRS(:,:,:)))
    WRITE(3,*)STEP,DABS(DP-DPOLD)
    
    
        
    !    REWIND 10
    !
    !    WRITE(10,*)"VARIABLES= X,Y,Z,DXY,U,V,W,SHEAR,RELAXATION"
    !    WRITE(10,*)"ZONE I=" ,XMAX , "J=",YMAX , "K=",ZMAX , "F=POINT"
    !
    !  
    !   
    !! ======== RHO CONTURE WRITING========
    !DO Z=1,ZMAX
    !  DO Y=1,YMAX
    !    DO X=1,XMAX
    !       WRITE(10,*)X,Y,Z,DXY(X,Y,Z),UPP(X,Y,Z),VPP(X,Y,Z),WPP(X,Y,Z),SHEAR(X,Y,Z),RELAXATION(X,Y,Z)
    !      END DO 
    !    END DO
    !END DO
    
    
    
    
   !   WRITE(FILE_NAME, "('output_',I5.5,'.dat')") STEP
   !   OPEN(UNIT=200, FILE=TRIM(FILE_NAME))
   !
   !   WRITE(200,*)"VARIABLES= X,Y,Z,Q_crit,U,V,W,Velocity-mag"
   !   WRITE(200,*)"ZONE I=" ,XMAX , "J=",YMAX , "K=",ZMAX , "F=POINT"
   !
   ! 
   !  
   ! ======== RHO CONTURE WRITING========
   !O Z=1,ZMAX
   ! DO Y=1,YMAX
   !   DO X=1,XMAX
   !     WRITE(200,*)X,Y,Z,Q_crit(X,Y,Z),UPP(X,Y,Z),VPP(X,Y,Z),WPP(X,Y,Z),USQRS(X,Y,Z)
   !     END DO 
   !   END DO
   !ND DO
   !
   !   CLOSE(UNIT=200)
     
    
    
    
    
    
    
    !!====== INDIVIDUAL VTK FILE EACH STEP FOR 3D ======
    !WRITE(MYFILE, "(a,i4.4,a)") "Output_", STEP, ".vtk"  ! Dynamically include timestep in the file name
    !OPEN(1, FILE=MYFILE, STATUS='NEW')  ! Ensure a new file is created for each timestep
    !
    !! Write VTK headers and dimensions
    !WRITE(1, "(A)") "# vtk DataFile Version 3.0"
    !WRITE(1, "(A)") "fluid_state"
    !WRITE(1, "(A)") "ASCII"
    !WRITE(1, "(A)") "DATASET RECTILINEAR_GRID"
    !WRITE(1, "(A,I5,I5,I5)") "DIMENSIONS", XMAX, YMAX, ZMAX
    !
    !! X_COORDINATES
    !WRITE(1, "(A,I5,A)") "X_COORDINATES ", XMAX, " float"
    !DO X = 0, XMAX - 1
    !    WRITE(1, "(I5)", advance='no') X + 1
    !END DO
    !WRITE(1, "(A)")  ! End line
    !
    !! Y_COORDINATES
    !WRITE(1, "(A,I5,A)") "Y_COORDINATES ", YMAX, " float"
    !DO Y = 0, YMAX - 1
    !    WRITE(1, "(I5)", advance='no') Y + 1
    !END DO
    !WRITE(1, "(A)")  ! End line
    !
    !! Z_COORDINATES
    !WRITE(1, "(A,I5,A)") "Z_COORDINATES ", ZMAX, " float"
    !DO Z = 0, ZMAX - 1
    !    WRITE(1, "(I5)", advance='no') Z + 1
    !END DO
    !WRITE(1, "(A)")  ! End line
    !
    !! POINT_DATA
    !WRITE(1, "(A,I6)") "POINT_DATA ", XMAX * YMAX * ZMAX
    !
    !! SCALARS Density float 1
    !WRITE(1, "(A)") "SCALARS Density float 1"
    !WRITE(1, "(A)") "LOOKUP_TABLE default"
    !DO Z = 1, ZMAX
    !    DO Y = 1, YMAX
    !        DO X = 1, XMAX
    !            WRITE(1, *) DXY(X, Y, Z)  ! 3D array for density
    !        END DO
    !    END DO
    !END DO
    !
    !! VECTORS Velocity float
    !WRITE(1, "(A)") "VECTORS Velocity float"
    !DO Z = 1, ZMAX
    !    DO Y = 1, YMAX
    !        DO X = 1, XMAX
    !            WRITE(1, *) UPP(X, Y, Z), VPP(X, Y, Z), WPP(X, Y, Z)  ! Add WPP for 3D velocity
    !        END DO
    !    END DO
    !END DO
    !
    !
    !
    !
    !
    !CLOSE(1)
    
    
    
    
    
    
    
    
    
    END SUBROUTINE