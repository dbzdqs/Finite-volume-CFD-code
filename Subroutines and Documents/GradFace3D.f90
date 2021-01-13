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
 Subroutine GradFace3D(Dim,Neq,NF,NP,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,Func,PrimFunc,DFuncX,DFuncY,DFuncZ)
 Implicit None
!*********************************************************************************************

 Integer::Dim,Neq,J,I,NP,NF,ME,NE,P1,P2,P3,P4,FacTyp,NFt
 Real(8)::vol_v,NXX,NYY,NZZ
 
 Real(8),Dimension(1:Neq,1:Dim)::Func,PrimFunc
 Real(8),Dimension(1:Neq,1:Dim)::DFuncX,DFuncY,DFuncZ
 
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:Dim)::XC,YC,ZC

 Integer,Dimension(1:4,1:8)::IDSt
 Integer,Dimension(1:8)::FaceTypet
 Real(8),Dimension(1:8)::NXt,NYt,NZt
 Real(8),Dimension(1:Neq)::Fu
!*********************************************************************************************

 DO I=1,NF
    FacTyp = FaceType(I)
 
    DFuncX(:,I) = 0.0
    DFuncY(:,I) = 0.0
    DFuncZ(:,I) = 0.0

    ME = IDS(1,I)
    NE = IDS(2,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
    P3 = IDS(5,I)
    P4 = IDS(6,I)

    X(NP+1) = XC(ME) ; Y(NP+1) = YC(ME) ; Z(NP+1) = ZC(ME)

    IF(NE/=0)THEN
     X(NP+2) = XC(NE) ; Y(NP+2) = YC(NE) ; Z(NP+2) = ZC(NE) 
    END IF

    FaceTypet(:)=3
    NFt = 0
    IF(FacTyp==3)Then
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P2 ; IDSt(3,NFt) = P1
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P3 ; IDSt(3,NFt) = P2
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P1 ; IDSt(3,NFt) = P3
    ElseIF(FacTyp==4)Then
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P2 ; IDSt(3,NFt) = P1
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P3 ; IDSt(3,NFt) = P2
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P4 ; IDSt(3,NFt) = P3
     NFt = NFt+1 ; IDSt(1,NFt) = NP+1 ; IDSt(2,NFt) = P1 ; IDSt(3,NFt) = P4
    EndIF

    IF(NE==0 .and. FacTyp==3)Then
     NFt = NFt+1 ; IDSt(1,NFt) = P1 ; IDSt(2,NFt) = P2 ; IDSt(3,NFt) = P3
    ElseIF(NE==0 .and. FacTyp==4)Then
     NFt = NFt+1 ; IDSt(1,NFt) = P1 ; IDSt(2,NFt) = P2 ; IDSt(3,NFt) = P3 ; IDSt(4,NFt) = P4
     FaceTypet(NFt)=4
    EndIF

    IF(NE/=0 .and. FacTyp==3)Then
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P1 ; IDSt(3,NFt) = P2
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P2 ; IDSt(3,NFt) = P3
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P3 ; IDSt(3,NFt) = P1
    ElseIF(NE/=0 .and. FacTyp==4)Then
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P1 ; IDSt(3,NFt) = P2
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P2 ; IDSt(3,NFt) = P3
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P3 ; IDSt(3,NFt) = P4
     NFt = NFt+1 ; IDSt(1,NFt) = NP+2 ; IDSt(2,NFt) = P4 ; IDSt(3,NFt) = P1
    EndIF

    Call GeoCalAnyShape3D(Dim,NFt,IDSt,X,Y,Z,FaceTypet,vol_v,Nxt,Nyt,Nzt)   

    Func(:,NP+1)= PrimFunc(:,ME)
    IF(NE/=0) Func(:,NP+2)= PrimFunc(:,NE)

    Do J=1,NFt

       NXX = NXt(J)/vol_v
       NYY = NYt(J)/vol_v
       NZZ = NZt(J)/vol_v
    
       P1 = IDSt(1,J)
       P2 = IDSt(2,J)
       P3 = IDSt(3,J)
       P4 = IDSt(4,J)

       IF(FaceTypet(J)==3)Then
        Fu(:) = ( Func(:,P1)+Func(:,P2)+Func(:,P3) )/3.0
       ElseIF(FaceTypet(J)==4)Then
        Fu(:) = ( Func(:,P1)+Func(:,P2)+Func(:,P3)+Func(:,P4) )/4.0
       EndIF

       DFuncX(:,I) = DFuncX(:,I) + Fu(:)*NXX 
       DFuncY(:,I) = DFuncY(:,I) + Fu(:)*NYY 
       DFuncZ(:,I) = DFuncZ(:,I) + Fu(:)*NZZ

    END DO

END DO  

!*********************************************************************************************
 End
!###########################################################################################