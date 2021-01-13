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
!// Date: May., 15, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KFi_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,WNP1,WTNP1,WB,WTB,DKX,DKY,DKZ,DFiX,DFiY,DFiZ)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,WNP1,WTNP1,WB,WTB
 Intent(Out  )::DKX,DKY,DKZ,DFiX,DFiY,DFiZ

 Integer::Dim,I,J,NP,NC,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,ME,NE,Pt,Neq
 Real(8)::K,Fi
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::X,Y,Z,XC,YC,ZC,DW,KP,FiP,DKX,DKY,DKZ,DFiX,DFiY,DFiZ
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType,NEC
 Real(8),Dimension(:,:),Allocatable ::Func,PrimFunc,DFuncX,DFuncY,DFuncZ
!************************************************************************************************
!Part 1:
 DO I=1,NP
    KP(I)  = 0.0
    FiP(I) = 0.0
    NEC(I) = 0
 END DO

Do I=NF1+1,NF2

    ME = IDS(1,I)
    NE = IDS(2,I)

    K  = 0.5 * ( WTNP1(1,ME)/WnP1(1,ME) + WTNP1(1,NE)/WnP1(1,NE) )
    Fi = 0.5 * ( WTNP1(2,ME)/WnP1(1,ME) + WTNP1(2,NE)/WnP1(1,NE) )
 
    DO J=1,FaceType(I)
       Pt=IDS(J+2,I) 
       KP(Pt)  = KP(Pt)  + K
       FiP(Pt) = FiP(Pt) + Fi
       NEC(Pt) = NEC(Pt) + 1
    End Do
 
End Do 

Do I=NF2+1,NF

    K  = WTB(1,I)/WB(1,I)
    Fi = WTB(2,I)/WB(1,I)
    
    DO J=1,FaceType(I)
        Pt=IDS(J+2,I) 
        KP(Pt)  = KP(Pt)  + K
        FiP(Pt) = FiP(Pt) + Fi
        NEC(Pt) = NEC(Pt) + 1
    End Do
 
End Do

 DO I=1,NP
    IF(NEC(I)==0)Cycle
    KP(I) = KP(I)/NEC(I)
    FiP(I) = FiP(I)/NEC(I)
 End Do

 DO I=NFW1+1,NFW2
     
    Fi = WTB(2,I)/WB(1,I)
    
    DO J=1,FaceType(I)
        Pt=IDS(J+2,I) 
        KP(Pt)  = 0.0
        FiP(Pt) = Fi
    End Do
 
 End Do

 
 Neq = 2
 Allocate( Func(1:Neq,1:Dim),PrimFunc(1:Neq,1:Dim) )
 Allocate( DFuncX(1:Neq,1:Dim),DFuncY(1:Neq,1:Dim),DFuncZ(1:Neq,1:Dim) )

 DO I=1,NP
    Func(1,I) = KP(I) 
    Func(2,I) = FiP(I) 
 End Do

 DO I=1,NC
    PrimFunc(1,I) = WTNP1(1,I)/WNP1(1,I)
    PrimFunc(2,I) = WTNP1(2,I)/WNP1(1,I)
 End Do

 Call GradFace3D(Dim,Neq,NF,NP,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,Func,PrimFunc,DFuncX,DFuncY,DFuncZ)

 Do I=1,NF
    DKX(I)  = DFuncX(1,I) ; DKY(I)  = DFuncY(1,I) ; DKZ(I)  = DFuncZ(1,I)
    DFiX(I) = DFuncX(2,I) ; DFiY(I) = DFuncY(2,I) ; DFiZ(I) = DFuncZ(2,I)
 End Do

  Do I=NFS1+1,NFS2
    DKZ(I)  = 0.0
    DFiZ(I) = 0.0
  End Do
  
 Deallocate( Func,PrimFunc,DFuncX,DFuncY,DFuncZ )
!*********************************************************************************************
 End
!###########################################################################################

