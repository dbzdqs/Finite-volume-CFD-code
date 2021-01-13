!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//       /////////////       ////////////////    ////////     //////    ////////////////  //!
!//       /////////////      ////////////////    //////////   //////    ////////////////   //!
!//      /////    /////     //////    //////    //////////// //////    /////               //!
!//     /////    //////    ////////////////    ///////////////////    ////////////////     //!
!//    /////    //////    ////////////////    ////// ////////////               /////      //!
!//   ///////////////    //////    //////    //////   //////////    ////////////////       //!
!// ///////////////     //////    //////    //////     ////////    ////////////////        //!
!//    Developer            Assistant    in      Numerical             Sciences            //!
!//----------------------------------------------------------------------------------------//!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Oct., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine VelTemp_GradFace3D(Dim,NC,NF,NP,NF1,NF2,NFW1,NFW2,NFS1,NFS2,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,&
                              WNP1,WB,GM,P,NX,NY,NZ,&
                              DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ)
 Implicit None
 !*********************************************************************************************
 Integer::Dim,J,I,k,NP,NF,NC,NF1,NF2,NFW1,NFW2,NFS1,NFS2,ME,NE,Pt,P1,P2,P3,P4,Neq
 Real(8)::GM,U,V,W,Temp ,Xm,Ym,ux,uy,vx,vy,NXX,NYY,NZZ
 Real(8),parameter::eps=0.00000001
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:Dim)::NEC
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::X,Y,Z,XC,YC,ZC,T,P
 Real(8),Dimension(1:Dim)::UP,VP,WP,TP,NX,NY,NZ
 Real(8),Dimension(1:Dim)::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ

 Real(8),Dimension(:,:),Allocatable ::Func,PrimFunc,DFuncX,DFuncY,DFuncZ
!*********************************************************************************************
!Part 1:
 Do J=1,NC
    T(J) = GM*P(J)/WNP1(1,J)         
 End Do

!Part 2:
 DO I=1,NP
    UP(I) = 0.0
    VP(I) = 0.0
    WP(I) = 0.0
    TP(I) = 0.0
    NEC(I)= 0
 END DO
 
!Part 3: 
 Do I=NF1+1,NF2
 
   !Part 4:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 5:
    U    = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )
    V    = 0.5*( WNP1(3,ME)/WNP1(1,ME) + WNP1(3,NE)/WNP1(1,NE) )
    W    = 0.5*( WNP1(4,ME)/WNP1(1,ME) + WNP1(4,NE)/WNP1(1,NE) )
    
    Temp = 0.5*( T(ME) + T(NE) )
    
   !Part 6:
    DO K=1,FaceType(I)
 
        Pt=IDS(K+2,I) 
 
        UP(Pt) = UP(Pt) + U
        VP(Pt) = VP(Pt) + V
        WP(Pt) = WP(Pt) + W    
        TP(Pt) = TP(Pt) + Temp
        NEC(Pt) = NEC(Pt) + 1

    End Do
 
 End Do 

!Part 7:
 Do I=NF2+1,NF

    U    = WB(2,I) /WB(1,I)
    V    = WB(3,I) /WB(1,I)
    W    = WB(4,I) /WB(1,I)
    Temp = GM*WB(6,I) /WB(1,I)
 
   !Part 8:
    DO K=1,FaceType(I)
  
        Pt=IDS(K+2,I) 
 
        UP(Pt) = UP(Pt) + U
        VP(Pt) = VP(Pt) + V
        WP(Pt) = WP(Pt) + W    
        TP(Pt) = TP(Pt) + Temp
        NEC(Pt) = NEC(Pt) + 1

    End Do
 
End Do

!Part 9:
 DO I=1,NP
    IF(NEC(I)==0)Cycle
    UP(I) = UP(I)/NEC(I)
    VP(I) = VP(I)/NEC(I)
    WP(I) = WP(I)/NEC(I)    
    TP(I) = TP(I)/NEC(I)
 END DO

!Part 10:
 DO I=NFW1+1,NFW2
    DO K=1,FaceType(I)
  
        Pt=IDS(K+2,I) 
    
        UP(Pt) = 0.0
        VP(Pt) = 0.0
        WP(Pt) = 0.0    
    
    END DO  
 END DO

!Part 11
 Neq = 4
  
!Part 12:
 Allocate( Func(1:Neq,1:Dim),PrimFunc(1:Neq,1:Dim) )
 Allocate( DFuncX(1:Neq,1:Dim),DFuncY(1:Neq,1:Dim),DFuncZ(1:Neq,1:Dim) )
 
!Part 13:
 DO I=1,NP
    Func(1,I) = UP(I) 
    Func(2,I) = VP(I) 
    Func(3,I) = WP(I) 
    Func(4,I) = TP(I)
 End Do

!Part 14:
 DO I=1,NC
    PrimFunc(1,I) = WNP1(2,I)/WNP1(1,I)
    PrimFunc(2,I) = WNP1(3,I)/WNP1(1,I)
    PrimFunc(3,I) = WNP1(4,I)/WNP1(1,I) 
    PrimFunc(4,I) = T(I)
 End Do

!Part 15:
 Call GradFace3D(Dim,Neq,NF,NP,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,Func,PrimFunc,DFuncX,DFuncY,DFuncZ)
      
!Part 16:
 Do I=1,NF
    DUX(I) = DFuncX(1,I) ; DUY(I) = DFuncY(1,I) ; DUZ(I) = DFuncZ(1,I)
    DVX(I) = DFuncX(2,I) ; DVY(I) = DFuncY(2,I) ; DVZ(I) = DFuncZ(2,I)
    DWX(I) = DFuncX(3,I) ; DWY(I) = DFuncY(3,I) ; DWZ(I) = DFuncZ(3,I)
    DTX(I) = DFuncX(4,I) ; DTY(I) = DFuncY(4,I) ; DTZ(I) = DFuncZ(4,I)
 End Do

!Part 17:
 Do I=NFS1+1,NFS2
  
    NXX = dabs(NX(I)) 
    NYY = dabs(NY(I))  
    NZZ = dabs(NZ(I))
    
    IF (NXX<eps .AND. NYY<eps) THEN
    DUZ(I) = 0.0
    DVZ(I) = 0.0
    DWZ(I) = 0.0
    DTZ(I) = 0.0
    END IF
    IF (NYY<eps .AND. NZZ<eps) THEN
    DUX(I) = 0.0
    DVX(I) = 0.0
    DWX(I) = 0.0
    DTX(I) = 0.0
    END IF
    IF (NXX<eps .AND. NZZ<eps) THEN
    DUY(I) = 0.0
    DVY(I) = 0.0
    DWY(I) = 0.0
    DTY(I) = 0.0
    END IF

 End Do
 
!Part 18: 
 Deallocate( Func,PrimFunc,DFuncX,DFuncY,DFuncZ )
!*********************************************************************************************
 End
!###########################################################################################
 

    
    
    
    
    
    
    
    

 
 
 
    
    
    
    
    
    
