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
!// Developed by: R. Amery, Mathmatical, Amirkabir university of Technology                //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine MoveMeshForBoundarylayer(Dim,NC,NR,NF,NP,NPtCurvs,NBL,NFR,BC,IDS,BLPt,XBL,YBL,X,Y)
 Implicit None
!************************************************************************************* 
 Intent(In   )::Dim,NR,NF,NP,NPtCurvs,NBL,NFR,BC,IDS,BLPt,XBL,YBL
 Intent(InOut)::X,Y,NC

 Integer::Dim,I,J,J1,JJ,II,P1,P2,P,NC
 Integer::NR
 Integer::NF !Number of Faces Constructing Mesh
 Integer::NP
 Integer::NPtCurvs
 Integer::NBL
 Integer::NBP
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:100)::BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
 Integer,Dimension(1:25,1:Dim)::BLPt
 Integer,Dimension(1:Dim)::IBP
 Real(8),Dimension(Dim)::DelX,DelY
 Real(8),Dimension(1:Dim)::XBL,YBL
 
 Integer,Dimension(1:Dim)::NConectPoints
 Integer,Dimension(1:10,1:Dim)::IConectPoints
!************************************************************************************* 
!Part 1:
 Call BoundPointLabeling(Dim,IDS,NR,NFR,BC,NBP,IBP)
 
 !!!do j=1,NBL
 !!!
 !!!   !Part 2:
 !!!    Do I=1,NPtCurvs
 !!!       P1 = BLPt(j    ,I)
 !!!       P2 = BLPt(j+1,I)
 !!!   
 !!!       P  = BLPt(NBL+2,I)
 !!! 
 !!!       DelX(P) = XBL(P2) - XBL(P1)
 !!!       DelY(P) = YBL(P2) - YBL(P1)
 !!!    End Do	
 !!!
 !!!   !Part 3:
 !!!    Call RBF_Moving_Mesh(Dim,NBP,NP,IBP,X,Y,DelX,DelY)
 !!!
 !!!   !Part 4:
 !!!    Do I=1,NP
 !!!       X(I) = X(I) + DelX(I)
 !!!       Y(I) = Y(I) + DelY(I)
 !!!    End Do
 !!!enddo
 
 DelX(:) = 0.0
 DelY(:) = 0.0

!Part 2:
 Do I=1,NPtCurvs
    P1 = BLPt(1    ,I)
    P2 = BLPt(NBL+1,I)
    
    P  = BLPt(NBL+2,I)
  
    DelX(P) = XBL(P2) - XBL(P1)
    DelY(P) = YBL(P2) - YBL(P1)
 End Do	
 
!Part 3:
 !!!Call RBF_Moving_Mesh(Dim,NBP,NP,IBP,X,Y,DelX,DelY)
 !call ConectPoints(Dim,NF,IDS,NConectPoints,IConectPoints)
 !call Linear_Spring2D(Dim,NConectPoints,IConectPoints,NP,X,Y,DelX,DelY)

!Part 4:
 Do I=1,NP
    X(I) = X(I) + DelX(I)
    Y(I) = Y(I) + DelY(I)
 End Do
!************************************************************************************  
 End
!####################################################################################