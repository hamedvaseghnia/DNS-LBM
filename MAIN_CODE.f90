!====================================================================================================================
!====================================================================================================================
!===      This Fortran Code is written for                                                                    =======
!===      Non-Newtonian Carreau-Yasuda fluid flow in Lattice Boltzmann Method                                 =======
!===      by Hamed Vaseghnia PhD candidate in computational engineering                                       =======
!===      UNIVERSITY of Stavanger                                                                             =======
!===      2023                                                                                                =======
!====================================================================================================================
!====================================================================================================================
!====================================================================================================================


PROGRAM SCMP
    USE VARIABLES
    IMPLICIT NONE
    
    
    NAMELIST/GEOMETRY/ XMAX,YMAX,ZMAX,RADIUS
    NAMELIST/PHYSICS/ DXYL,DXYV
    NAMELIST/MULTIPHASE_METHOD/S_C,R_K,P_R,C_S
    NAMELIST/SC/GE,DXYC
    NAMELIST/RK/A,B,R,TR,TC,DXYC
    NAMELIST/PR/A,B,R,AF,TR,TC,DXYC,KAPPA,SIGMA
    NAMELIST/CSS/A,B,R,TR,TC,DXYC
    
    OPEN(100,FILE="PARAMETERS.PAR")
    READ(100,GEOMETRY)
    READ(100,PHYSICS)
    READ(100,MULTIPHASE_METHOD)
    
    !CALL SYSTEM_CLOCK(count_rate=rate)
    
    PRINT*," ITTERATED:","  INSIDE PRESSURE-OUTSIDE PRESSURE="
    
        ALLOCATE(U(XMAX,YMAX,ZMAX))
        ALLOCATE(V(XMAX,YMAX,ZMAX))
        ALLOCATE(W(XMAX,YMAX,ZMAX))
        ALLOCATE(DXY(XMAX,YMAX,ZMAX))
        ALLOCATE(UXY(XMAX,YMAX,ZMAX))
        ALLOCATE(USQR(XMAX,YMAX,ZMAX))
        ALLOCATE(VU(XMAX,YMAX,ZMAX))
        ALLOCATE(VV(XMAX,YMAX,ZMAX))
        ALLOCATE(VOR(XMAX,YMAX,ZMAX))
        ALLOCATE(FIN(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(FEQ(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(FOUT(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(FX_GOU(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(FY_GOU(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(FZ_GOU(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(S(XMAX,YMAX,ZMAX))
        ALLOCATE(F_SCX(XMAX,YMAX,ZMAX))
        ALLOCATE(F_SCY(XMAX,YMAX,ZMAX))
        ALLOCATE(F_SCZ(XMAX,YMAX,ZMAX))
        ALLOCATE(P(XMAX,YMAX,ZMAX))
        ALLOCATE(WALL(XMAX,YMAX,ZMAX))
        ALLOCATE(UP(XMAX,YMAX,ZMAX))
        ALLOCATE(VP(XMAX,YMAX,ZMAX))
        ALLOCATE(WP(XMAX,YMAX,ZMAX))
        ALLOCATE(VPP(XMAX,YMAX,ZMAX))
        ALLOCATE(UPP(XMAX,YMAX,ZMAX))
        ALLOCATE(WPP(XMAX,YMAX,ZMAX))
        ALLOCATE(USQRS(XMAX,YMAX,ZMAX))
        !ALLOCATE(GIN(0:6,XMAX,YMAX,ZMAX))
        !ALLOCATE(GEQ(0:6,XMAX,YMAX,ZMAX))
        !ALLOCATE(GOUT(0:6,XMAX,YMAX,ZMAX))
        ALLOCATE(TXY(XMAX,YMAX,ZMAX))
        ALLOCATE(UREAL(XMAX,YMAX,ZMAX))
        ALLOCATE(VREAL(XMAX,YMAX,ZMAX))
        ALLOCATE(WREAL(XMAX,YMAX,ZMAX))
        !ALLOCATE(DPDT(XMAX,YMAX,ZMAX))
        !ALLOCATE(DIVERGENCE(XMAX,YMAX,ZMAX))
        !ALLOCATE(TGRAD(XMAX,YMAX,ZMAX))
        !ALLOCATE(TDIVER(XMAX,YMAX,ZMAX))
        ALLOCATE(MFIN(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(MFEQ(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(TRANSFORMED(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(RELAXEDM(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(NEQ(0:18,XMAX,YMAX,ZMAX))
        !ALLOCATE(MGIN(0:6,XMAX,YMAX,ZMAX))
        !ALLOCATE(MGEQ(0:6,XMAX,YMAX,ZMAX))
        !ALLOCATE(RELAXEDN(0:6,XMAX,YMAX,ZMAX))
        !ALLOCATE(TRANSFORMEDN(0:6,XMAX,YMAX,ZMAX))
        ALLOCATE(TMP_RST(XMAX,YMAX,ZMAX))
        ALLOCATE(MFX_GOU(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(MFY_GOU(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(MFZ_GOU(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(RELAXEDM_FX(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(RELAXEDM_FY(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(RELAXEDM_FZ(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(MRT_FX(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(MRT_FY(0:18,XMAX,YMAX,ZMAX))
        ALLOCATE(MRT_FZ(0:18,XMAX,YMAX,ZMAX))
        
        ALLOCATE(Sab_XX(XMAX,YMAX,ZMAX))
        ALLOCATE(Sab_YY(XMAX,YMAX,ZMAX))
        ALLOCATE(Sab_ZZ(XMAX,YMAX,ZMAX))
        ALLOCATE(Sab_XY(XMAX,YMAX,ZMAX))
        ALLOCATE(Sab_XZ(XMAX,YMAX,ZMAX))
        ALLOCATE(Sab_YZ(XMAX,YMAX,ZMAX))
        
        ALLOCATE(RST_XX(XMAX,YMAX,ZMAX))
        ALLOCATE(RST_YY(XMAX,YMAX,ZMAX))
        ALLOCATE(RST_XY(XMAX,YMAX,ZMAX))
        ALLOCATE(RST_XZ(XMAX,YMAX,ZMAX))
        ALLOCATE(RST_YZ(XMAX,YMAX,ZMAX))
        ALLOCATE(RST_ZZ(XMAX,YMAX,ZMAX))
        
        ALLOCATE(INVARIANT(XMAX,YMAX,ZMAX))
        ALLOCATE(SHEAR(XMAX,YMAX,ZMAX))
        ALLOCATE(RELAXATION(XMAX,YMAX,ZMAX))
        ALLOCATE(VISCOSITY(XMAX,YMAX,ZMAX))
        
        ALLOCATE(SS(0:18,XMAX,YMAX,ZMAX))
        
        ALLOCATE(S_11(XMAX,YMAX,ZMAX))
        ALLOCATE(S_12(XMAX,YMAX,ZMAX))
        ALLOCATE(S_13(XMAX,YMAX,ZMAX))
        ALLOCATE(S_22(XMAX,YMAX,ZMAX))
        ALLOCATE(S_23(XMAX,YMAX,ZMAX))
        ALLOCATE(S_33(XMAX,YMAX,ZMAX))
        
        ALLOCATE(S_SQUARED(XMAX,YMAX,ZMAX))
        ALLOCATE(omega_SQUARED(XMAX,YMAX,ZMAX))
        ALLOCATE(Q_crit(XMAX,YMAX,ZMAX))
        ALLOCATE(local_energy(XMAX,YMAX,ZMAX))
        
        
    
    
        
        
        
    OPEN(10,FILE="Kinetic energy Normalized.txt")
    OPEN(11,FILE="Kinetic energy.txt")
    !OPEN(11,FILE="VELOCITY PROFILE.csv")
    !OPEN(40,FILE="CPU RUNTIME.TXT")
    !!OPEN(200,FILE="RHO OUTPUT_Physical.TEC")
    !OPEN(400,FILE="WSS.CSV")
    
    
    CALL WALL_COORDINATES
    
    CALL INITIALIZE
    
    CALL INITDENSITY
    
    print*,"initializtion done"
    energy_diss_old=0d0
    
    DPOLD=0.D0
    DO STEP=1,70000  !========= INITIAL LOOP STARTS =========== 
    !CALL SYSTEM_CLOCK(start)
    
    CALL MACROS
    
    CALL F_SC
    
    CALL F_EQ 
    
    CALL COLLISION
    
    CALL STREAMING
    
    !CALL BOUNCEBACK
    
    !CALL PRESS_BOUND
    
    DO Z=1,ZMAX
        DO Y=1,YMAX
            DO X=1,XMAX
    
                UPP(X,Y,Z)=U(X,Y,Z)+0.5D0*F_SCX(X,Y,Z)/DXY(X,Y,Z)
                VPP(X,Y,Z)=V(X,Y,Z)+0.5D0*F_SCY(X,Y,Z)/DXY(X,Y,Z)
                WPP(X,Y,Z)=W(X,Y,Z)+0.5D0*F_SCZ(X,Y,Z)/DXY(X,Y,Z)
                USQRS(X,Y,Z)=dsqrt(UPP(X,Y,Z)*UPP(X,Y,Z)+VPP(X,Y,Z)*VPP(X,Y,Z)+WPP(X,Y,Z)*WPP(X,Y,Z))
    
            END DO 
        END DO
    END DO
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    energy_diss = 0d0
    
    DO Z = 1, ZMAX
        DO Y = 1, YMAX
            DO X = 1, XMAX
                
    
                    local_energy(X,Y,Z)= 0.5d0*dxy(x,y,z)*(UREAL(X,Y,Z)**2d0 + VREAL(X,Y,Z)**2d0 + WREAL(X,Y,Z)**2d0)
                
            END DO
        END DO
    END DO
    
    
    
    DO Z = 1, ZMAX
        DO Y = 1, YMAX
            DO X = 1, XMAX
                
                energy_diss = energy_diss + local_energy(X,Y,Z)
                
            END DO
        END DO
    END DO
    
    energy_diss=energy_diss/(xmax*ymax*zmax)
    
    !write(10,*) 2d0*pi*0.1018591d0*step/ymax, energy_diss/(0.125*0.1018591d0**2)
    !write(10,*) 2d0*pi*0.1018591d0*step/ymax, energy_diss
    
    dissipation_rate=-(energy_diss - energy_diss_old)
    energy_diss_old=energy_diss
    write(11,*)step,energy_diss
    write(10,*)step,dissipation_rate
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
     IF (MOD(STEP,100)==0) THEN
    
     CALL OUTPUT
     
      
     PRINT*, STEP,DABS(DP-DPOLD)
     IF (DABS(DP-DPOLD) .LT. 1.D-8)THEN
         
        
     PRINT*,"SOLUTION IS CONVERGED"
     EXIT
     ELSE
     DPOLD=DP
     END IF
     END IF
    
     
    END DO             !=========  LOOP ENDS ============
    
    
     PRINT *,"TIME=",(FINISH-START)/60,"MINUTES"
     PAUSE
    END PROGRAM SCMP