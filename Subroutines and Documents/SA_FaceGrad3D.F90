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
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine SA_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,NX,NY,NZ,&
                          WNP1,WTNP1,WB,WTB,DNuX,DNuY,DNuZ)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,WNP1,WTNP1,WB,WTB
 Intent(Out  )::DNuX,DNuY,DNuZ

 Integer::Dim,I,J,NP,NC,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,ME,NE,Pt,Neq
 Real(8)::NXX,NYY,NZZ,Nu
 Real(8),parameter::eps=0.00000001
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::X,Y,Z,XC,YC,ZC,NX,NY,NZ,DW,NuP,DNuX,DNuY,DNuZ
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType,NEC
 Real(8),Dimension(:,:),Allocatable ::Func,PrimFunc,DFuncX,DFuncY,DFuncZ
!************************************************************************************************
!Part 1:
 DO I=1,NP
    NuP(I) = 0.0
    NEC(I) = 0
 END DO

Do I=NF1+1,NF2

    ME = IDS(1,I)
    NE = IDS(2,I)

    Nu = 0.5 * ( WTNP1(1,ME)/WnP1(1,ME) + WTNP1(1,NE)/WnP1(1,NE) )
 
    DO J=1,FaceType(I)
       Pt=IDS(J+2,I) 
       NuP(Pt)  = NuP(Pt)  + Nu
       NEC(Pt) = NEC(Pt) + 1
    End Do
 
End Do 

Do I=NF2+1,NF

    Nu  = WTB(1,I)/WB(1,I)
    
    DO J=1,FaceType(I)
        Pt=IDS(J+2,I) 
        NuP(Pt) = NuP(Pt) + Nu
        NEC(Pt) = NEC(Pt) + 1
    End Do
 
End Do

 DO I=1,NP
    IF(NEC(I)==0)Cycle
    NuP(I)  = NuP(I)/NEC(I)
 End Do

 DO I=NFW1+1,NFW2
     
    Nu = WTB(1,I)/WB(1,I)
    
    DO J=1,FaceType(I)
        Pt=IDS(J+2,I) 
        NuP(Pt)  = 0.0
    End Do
 
 End Do

 
 Neq = 1
 Allocate( Func(1:Neq,1:Dim),PrimFunc(1:Neq,1:Dim) )
 Allocate( DFuncX(1:Neq,1:Dim),DFuncY(1:Neq,1:Dim),DFuncZ(1:Neq,1:Dim) )

 DO I=1,NP
    Func(1,I) = NuP(I) 
 End Do

 DO I=1,NC
    PrimFunc(1,I) = WTNP1(1,I)/WNP1(1,I)
 End Do

 Call GradFace3D(Dim,Neq,NF,NP,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,Func,PrimFunc,DFuncX,DFuncY,DFuncZ)

 Do I=1,NF
    DNuX(I)  = DFuncX(1,I) ; DNuY(I)  = DFuncY(1,I) ; DNuZ(I)  = DFuncZ(1,I)
 End Do

  Do I=NFS1+1,NFS2
  
    NXX = dabs(NX(I)) 
    NYY = dabs(NY(I)) 
    NZZ = dabs(NZ(I))
    
    IF (NXX<eps .AND. NYY<eps) THEN
    DNuZ(I)  = 0.0
    END IF
    IF (NYY<eps .AND. NZZ<eps) THEN
    DNuX(I)  = 0.0
    END IF
    IF (NXX<eps .AND. NZZ<eps) THEN
    DNuY(I)  = 0.0
    END IF

  End Do
  
 Deallocate( Func,PrimFunc,DFuncX,DFuncY,DFuncZ )
!*********************************************************************************************
 End
!###########################################################################################

