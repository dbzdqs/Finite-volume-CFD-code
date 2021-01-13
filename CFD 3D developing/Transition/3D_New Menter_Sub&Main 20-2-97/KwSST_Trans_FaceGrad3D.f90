!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the gradient at faces of Mesh for 2 eq. Turbulence Model 2D   //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F054F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KwSST_Trans_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,NX,NY,NZ,&
                                   WNP1,WTNP1,WB,WTB,DKX,DKY,DKZ,DFiX,DFiY,DFiZ,DGamaX,DGamaY,DGamaZ)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,NX,NY,NZ,WNP1,WTNP1,WB,WTB
 Intent(Out  )::DKX,DKY,DKZ,DFiX,DFiY,DFiZ,DGamaX,DGamaY,DGamaZ

 Integer::Dim,I,J,NP,NC,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,ME,NE,Pt,Neq
 Real(8)::NXX,NYY,NZZ,K,Fi,Gam
 Real(8),parameter::eps=0.00000001
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::X,Y,Z,XC,YC,ZC,NX,NY,NZ,DW,KP,FiP,GP,DKX,DKY,DKZ,DFiX,DFiY,DFiZ,DGamaX,DGamaY,DGamaZ
 Real(8),Dimension(1:3,1:Dim)::WTNP1,WTB
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType,NEC
 Real(8),Dimension(:,:),Allocatable ::Func,PrimFunc,DFuncX,DFuncY,DFuncZ
!************************************************************************************************
!Part 1:
 DO I=1,NP
    KP(I)  = 0.0
    FiP(I) = 0.0
    GP(I) = 0.0
    NEC(I) = 0
 END DO

!PART 2:
 Do I=NF1+1,NF2

    !PART 3:
    ME = IDS(1,I)
    NE = IDS(2,I)
   
   !PART 4:
    K  = 0.5 *  ( WTNP1(1,ME)/WnP1(1,ME) + WTNP1(1,NE)/WnP1(1,NE) )
    Fi = 0.5 *  ( WTNP1(2,ME)/WnP1(1,ME) + WTNP1(2,NE)/WnP1(1,NE) )
    Gam = 0.5* ( WTNP1(3,ME)/WnP1(1,ME) + WTNP1(3,NE)/WnP1(1,NE) )
 
   !PART 5:
    DO J=1,FaceType(I)
       Pt=IDS(J+2,I) 
       KP(Pt)  = KP(Pt)  + K
       FiP(Pt) = FiP(Pt) + Fi
       GP(Pt)  = GP(Pt)  + Gam
       NEC(Pt) = NEC(Pt) + 1
    End Do
 
 End Do 
 
!PART 6:
Do I=NF2+1,NF

   !PART 7:
    K  = WTB(1,I)/WB(1,I)
    Fi = WTB(2,I)/WB(1,I)
    Gam = WTB(3,I)/WB(1,I)
    
   !PART 8:
    DO J=1,FaceType(I)
        Pt=IDS(J+2,I) 
        KP(Pt)  = KP(Pt)  + K
        FiP(Pt) = FiP(Pt) + Fi
        GP(Pt)  = GP(Pt) + Gam
        NEC(Pt) = NEC(Pt) + 1
    End Do
 
End Do

!PART 9:
 DO I=1,NP
    IF(NEC(I)==0)Cycle
    KP(I)  = KP(I)/NEC(I)
    FiP(I) = FiP(I)/NEC(I)
    GP(I)  = GP(I)/NEC(I)
 End Do
 
!PART 10:
 DO I=NFW1+1,NFW2
    
   !PART 11:
    Fi = WTB(2,I)/WB(1,I)
    Gam = WTB(3,I)/WB(1,I)
    
   !PART 12:
    DO J=1,FaceType(I)
        Pt=IDS(J+2,I) 
        KP(Pt)  = 0.0
        FiP(Pt) = Fi
        GP(Pt)  = Gam
    End Do
 
 End Do

!PART 13:
 Neq = 3
 
!PART 14:
 Allocate( Func(1:Neq,1:Dim),PrimFunc(1:Neq,1:Dim) )
 Allocate( DFuncX(1:Neq,1:Dim),DFuncY(1:Neq,1:Dim),DFuncZ(1:Neq,1:Dim) )

!PART 15:
 DO I=1,NP
    Func(1,I) = KP(I) 
    Func(2,I) = FiP(I)
    Func(3,I) = GP(I)
 End Do

!PART 16:
 DO I=1,NC
    PrimFunc(1,I) = WTNP1(1,I)/WNP1(1,I)
    PrimFunc(2,I) = WTNP1(2,I)/WNP1(1,I)
    PrimFunc(3,I) = WTNP1(3,I)/WNP1(1,I)
 End Do

!PART 17:
 Call GradFace3D(Dim,Neq,NF,NP,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,Func,PrimFunc,DFuncX,DFuncY,DFuncZ)

!PART 18:
 Do I=1,NF
    DKX(I)  = DFuncX(1,I) ; DKY(I)  = DFuncY(1,I) ; DKZ(I)  = DFuncZ(1,I)
    DFiX(I) = DFuncX(2,I) ; DFiY(I) = DFuncY(2,I) ; DFiZ(I) = DFuncZ(2,I)
    DGamaX(I) = DFuncX(3,I) ; DGamaY(I) = DFuncY(3,I) ; DGamaZ(I) = DFuncZ(3,I)
 End Do

 !PART 19:
  Do I=NFS1+1,NFS2
  
    NXX = dabs(NX(I)) 
    NYY = dabs(NY(I)) 
    NZZ = dabs(NZ(I))
    
    IF (NXX<eps .AND. NYY<eps) THEN
    DKZ(I)  = 0.0
    DFiZ(I) = 0.0
    DGamaZ(I) = 0.0
    END IF
    IF (NYY<eps .AND. NZZ<eps) THEN
    DKX(I)  = 0.0
    DFiX(I) = 0.0
    DGamaX(I) = 0.0
    END IF
    IF (NXX<eps .AND. NZZ<eps) THEN
    DKY(I)  = 0.0
    DFiY(I) = 0.0
    DGamaY(I) = 0.0
    END IF

  End Do
  
!PART 20:
 Deallocate( Func,PrimFunc,DFuncX,DFuncY,DFuncZ )
!*******************************************************************************************
 End
!###########################################################################################

