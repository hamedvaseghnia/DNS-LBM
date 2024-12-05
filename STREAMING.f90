SUBROUTINE STREAMING
    USE VARIABLES
    

    !!$OMP PARALLEL DO PRIVATE(X,Y,Z,XR,XL,YT,YB,ZI,ZO) SHARED(FIN,FOUT)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(X,Y,Z,XR,XL,YT,YB,ZI,ZO) SHARED(FIN,FOUT) SCHEDULE(STATIC)
    DO Z=1,ZMAX
        DO Y=1,YMAX
            DO X=1,XMAX
            
                YT = MOD(Y,YMAX) + 1 
                XR = MOD(X,XMAX) + 1 
                YB = YMAX - MOD(YMAX + 1 - Y, YMAX) 
                XL = XMAX - MOD(XMAX + 1 - X, XMAX)
                ZI = MOD(Z,ZMAX) + 1 
                ZO = ZMAX - MOD(ZMAX + 1 - Z, ZMAX)
    
                FIN(0,X,Y,Z)     = FOUT(0,X,Y,Z)
                FIN(1,XR,Y,Z)    = FOUT(1,X,Y,Z)
                FIN(2,XL,Y,Z)    = FOUT(2,X,Y,Z)
                FIN(3,X,YT,Z)    = FOUT(3,X,Y,Z)
                FIN(4,X,YB,Z)    = FOUT(4,X,Y,Z)
                FIN(5,X,Y,ZI)    = FOUT(5,X,Y,Z)
                FIN(6,X,Y,ZO)    = FOUT(6,X,Y,Z)
                FIN(7,XR,YT,Z)   = FOUT(7,X,Y,Z)
                FIN(8,XL,YB,Z)   = FOUT(8,X,Y,Z)
                FIN(9,XR,Y,ZI)   = FOUT(9,X,Y,Z)
                FIN(10,XL,Y,ZO)  = FOUT(10,X,Y,Z)
                FIN(11,X,YT,ZI)  = FOUT(11,X,Y,Z)
                FIN(12,X,YB,ZO)  = FOUT(12,X,Y,Z)
                FIN(13,XR,YB,Z)  = FOUT(13,X,Y,Z)
                FIN(14,XL,YT,Z)  = FOUT(14,X,Y,Z)
                FIN(15,XR,Y,ZO)  = FOUT(15,X,Y,Z)
                FIN(16,XL,Y,ZI)  = FOUT(16,X,Y,Z)
                FIN(17,X,YT,ZO)  = FOUT(17,X,Y,Z)
                FIN(18,X,YB,ZI)  = FOUT(18,X,Y,Z)
    
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    RETURN
    END SUBROUTINE
    