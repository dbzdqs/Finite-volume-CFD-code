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
 Subroutine RemoveRegn2D(Dim,LiveRgn,DeadRgn,NP_BL,NF_BL,NR_BL,NFR_BL,IDS_BL,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,LiveRgn,DeadRgn,NP_BL,X,Y
 Intent(InOut)::NF_BL,NR_BL,NFR_BL,IDS_BL

 Integer::Dim,I,J,NF,Heir,Dead,Live,Edge,P1P1,P2P2,P2P1,P1P2,P1_1,P2_1,P1_2,P2_2,LastFace,FirstFace
 Real(8)::Eps
 Integer::LiveRgn
 Integer::DeadRgn
 Integer::NP_BL
 Integer::NF_BL
 Integer::NR_BL
 Integer,Dimension(1:100)::NFR_BL
 Integer,Dimension(1:4,1:Dim)::IDS_BL
 Real(8),Dimension(1:Dim)::X,Y

 Integer,Dimension(1:Dim)::NConectEdge
 Integer,Dimension(1:Dim,1:100)::IConectEdge

 Integer::NdeadPt
 Integer,Dimension(1:Dim)::LivePt,DeadPt
 Integer::NdeadEdg
 Integer,Dimension(1:Dim)::LiveEdg,DeadEdg
!********************************************************************************************* 
!Part 1:
 eps=0.0001
 NdeadEdg = 0
 NdeadPt  = 0 
 
 FirstFace = 0
 Do I=1,LiveRgn-1
    FirstFace = FirstFace + NFR_BL(I)
 End Do
 LastFace = FirstFace + NFR_BL(LiveRgn)
 
   !print*,FirstFace+1,LastFace
   !print*,LastFace+1,LastFace+NFR_BL(LiveRgn) 
   ! pause
!Part 2:            
 Do i=FirstFace+1,LastFace
    P1_1 = IDS_BL(3,i) 
    P2_1 = IDS_BL(4,i)

   !Pa  rt 3:  
    Do j=LastFace+1,LastFace+NFR_BL(DeadRgn) 
       P1_2 = IDS_BL(3,j) 
       P2_2 = IDS_BL(4,j)
      
      !Part 4:  
       P1P1=0
       P2P2=0
       P1P2=0
       P2P1=0
          
       if( abs( x(P1_1)-x(P1_2) )<eps .and. abs( y(P1_1)-y(P1_2) )<eps ) P1P1=1
       if( abs( x(P2_1)-x(P2_2) )<eps .and. abs( y(P2_1)-y(P2_2) )<eps ) P2P2=1
      
       if( abs( x(P1_1)-x(P2_2) )<eps .and. abs( y(P1_1)-y(P2_2) )<eps ) P1P2=1
       if( abs( x(P2_1)-x(P1_2) )<eps .and. abs( y(P2_1)-y(P1_2) )<eps ) P2P1=1
       
    !print*,P1_1,P1_2,i,j  !abs( x(P1_1)-x(P1_2) ),abs( y(P1_1)-y(P1_2) ),P1P1
    !print*,P2_1,P2_2,i,j !abs( x(P2_1)-x(P2_2) ),abs( y(P2_1)-y(P2_2) ),P2P2
    !print*,P1_1,P2_2,i,j  !abs( x(P1_1)-x(P2_2) ),abs( y(P1_1)-y(P2_2) ),P1P2
    !print*,P2_1,P1_2,i,j  !abs( x(P2_1)-x(P1_2) ),abs( y(P2_1)-y(P1_2) ),P2P1
    !pause
      !Part 5:  
       if( P1P1==1 .and. P2P2==1 )Then
        NdeadEdg = NdeadEdg + 1
        LiveEdg(NdeadEdg) = i
        DeadEdg(NdeadEdg) = j
          
        NdeadPt =  NdeadPt + 1
        LivePt(NdeadPt) = P1_1
        DeadPt(NdeadPt) = P1_2

        NdeadPt =  NdeadPt + 1
        LivePt(NdeadPt) = P2_1
        DeadPt(NdeadPt) = P2_2

        exit
       endif
            
      !Part 6:     
       if( P2P1==1 .and. P1P2==1 )Then
        NdeadEdg = NdeadEdg + 1
        LiveEdg(NdeadEdg) = i
        DeadEdg(NdeadEdg) = j
         
        NdeadPt =  NdeadPt + 1
        LivePt(NdeadPt) = P1_1
        DeadPt(NdeadPt) = P2_2

        NdeadPt =  NdeadPt + 1
        LivePt(NdeadPt) = P2_1
        DeadPt(NdeadPt) = P1_2

        exit
       endif
      
      End Do 
  
 End Do 
 print*,'NdeadEdg',NdeadEdg
!Part 7:
 Do i=1,NdeadEdg  
    Live = LiveEdg(NdeadEdg)
    Dead = DeadEdg(NdeadEdg)
    IDS_BL(2,Live) = IDS_BL(1,Dead)
 End Do
 
!Part 8:
 Call ConectedEdgeOfPoint(Dim,NF_BL,NP_BL,IDS_BL,NConectEdge,IConectEdge)
 
!Part 9:  
 Do J=1,NDeadPt
    Live = LivePt(j)
    Dead = DeadPt(j)
     
     Do I=1,NConectEdge(Dead)
        Edge = IConectEdge(Dead,I)
        IF(IDS_BL(3,Edge)==Dead) IDS_BL(3,Edge) = Live
        IF(IDS_BL(4,Edge)==Dead) IDS_BL(4,Edge) = Live
     END Do
 End Do

!Part 10:
 NF_BL = NF_BL - NdeadEdg
 
!Part 11:
 Do i=NFR_BL(LiveRgn) + 1 , NF_BL   
    IDS_BL(:,i) = IDS_BL(:,i+NFR_BL(DeadRgn))
 End Do  

!Part 12:   
 Do i=DeadRgn,NR_BL   
    NFR_BL(i) = NFR_BL(i+1)
 End Do 
  
!Part 13:
 NR_BL = NR_BL - 1 
!*********************************************************************************************
 End
!###########################################################################################