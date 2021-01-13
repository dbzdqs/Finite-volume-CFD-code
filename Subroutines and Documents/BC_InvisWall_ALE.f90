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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BC_InvisWall_ALE(Dim,NFW1,NFW2,IDS,GM,P,NX,NY,DA,WNP1,Face_Velocity,WB)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,GM,P,NX,NY,DA,WNP1 ,Face_Velocity
 Intent(Out  )::WB

 Integer::Dim,I,NFW1,NFW2,ME
 Real(8)::GM,GM1,PB,U,V,V_Tang,TX,TY,NXX,NYY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::P,NX,NY ,DA ,GF
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:2,1:Dim)::Face_Velocity
!*********************************************************************************************	
!Part 1:
 GM1= GM-1.
	
!Part 2:
 DO I=NFW1+1,NFW2
	
   !Part 3:
    ME = IDS(1,I)
    NXX = NX(I)/DA(I)
    NYY = NY(I)/DA(I)
	
   !Part 4:
    WB(1,I) = WNP1(1,ME)
    
   !Part 5:
    TX = NYY
    TY =-NXX
	
   !Part 6:
    U = WNP1(2,ME) / WNP1(1,ME)
    V = WNP1(3,ME) / WNP1(1,ME)
    V_Tang = U*TX + V*TY
    
   !Part 7:
    WB(2,I) = WB(1,I) * ( V_Tang * TX + Face_Velocity(1,I) )
    WB(3,I) = WB(1,I) * ( V_Tang * TY + Face_Velocity(2,I) )
	
   !Part 8:
    WB(5,I) = P(ME)
    WB(4,I) =  WB(5,I)/GM1 + 0.5*( WB(2,I)* WB(2,I) + WB(3,I)*WB(3,I))/WB(1,I)

 End do
!*********************************************************************************************
 End
!###########################################################################################
