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
 Subroutine SA_BC(Dim,NFW1,NFW2,NF,IDS,NX,NY,WB,WTNP1,Nuhat0,WTB)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,NF,IDS,NX,NY,WB,WTNP1,Nuhat0
 Intent(Out  )::WTB

 Integer::Dim,I,P1,P2,ME,NFW1,NFW2,NF
 Real(8)::U,V,DX,DY,Un,Nuhat0
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,WTB,WTNP1
 Real(8),Dimension(1:5,1:Dim)::WB
!*********************************************************************************************
!Part 1:	
 Do I=NFW2+1,NF
	
   !Part 2:
    ME = IDS(1,I)

   !Part 3:
    U = WB(2,I)/WB(1,I) 
    V = WB(3,I)/WB(1,I) 
	
   !Part 4: 
    Un = U*NX(I) + V*NY(I)
	
   !Part 5:
	IF( Un<0. )then
	 WTB(I) = Nuhat0
    Else
	 WTB(I) = WTNP1(ME)
	Endif

 End do

!Part 6:
 DO I=NFW1+1,NFW2
    WTB(I) = 0.0
 End do	
!*********************************************************************************************
 End
!###########################################################################################
