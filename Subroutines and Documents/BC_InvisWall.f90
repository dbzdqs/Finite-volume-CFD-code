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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BC_InvisWall(Dim,NFW1,NFW2,IDS,GM,P,NX,NY,WNP1,WB)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,GM,P,NX,NY,WNP1
 Intent(Out  )::WB

 Integer::Dim,I,NFW1,NFW2,ME
 Real(8)::GM,GM1,PB,U,V,Temp,TX,TY,NXX,NYY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::P,NX,NY
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:4,1:Dim)::WNP1 
!*********************************************************************************************	
!Part 1:
 GM1= GM-1.
	
!Part 2:
 DO I=NFW1+1,NFW2
	
   !Part 3:
    ME = IDS(1,I)
	NXX = NX(I)
    NYY = NY(I)
	
   !Part 4:
    WB(1,I) = WNP1(1,ME)
	
   !Part 5:
    U = WNP1(2,ME) / WNP1(1,ME)
    V = WNP1(3,ME) / WNP1(1,ME)
    Temp = U*NYY - V*NXX
    
    IF (Temp>=0.0) Then
     TX = NYY
     TY =-NXX
    Else 
     TX =-NYY
     TY = NXX
    End IF
    	
   !Part 6:
    WB(2,I) = WB(1,I) * ( U*TX + V*TY ) * TX
    WB(3,I) = WB(1,I) * ( U*TX + V*TY ) * TY
	
   !Part 7:
	WB(5,I) = P(ME)
    WB(4,I) =  WB(5,I)/GM1 + 0.5*( WB(2,I)* WB(2,I) + WB(3,I)*WB(3,I))/WB(1,I)

 End do
!*********************************************************************************************
 End
!###########################################################################################
