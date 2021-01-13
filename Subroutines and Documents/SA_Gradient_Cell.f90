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
 Subroutine SA_Gradient_Cell(Dim,NC,NF1,NF2,NF,NX,NY,A,WTNP1,WNP1,WB,WTB,IDS,DNuXC,DNuYC,DRNuXC,DRNuYC)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,NX,NY,A,WTNP1,WNP1,WB,WTB,IDS
 Intent(Out  )::DNuXC,DNuYC,DRNuXC,DRNuYC

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,P1,P2
 Real(8)::NXX,NYY,AA,U,V,Nu,RNu,DNuXCj,DNuYCj
 Real(8),Dimension(1:Dim)::NX,NY,A,DNuXC,DNuYC,DRNuXC,DRNuYC,WTNP1,WTB
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:4,1:Dim)::IDS
!*********************************************************************************************
!Part 1:	
 DO I=1,NC
    DNuXC(I)  = 0.0
    DNuYC(I)  = 0.0

    DRNuXC(I) = 0.0
    DRNuYC(I) = 0.0

 End Do

!Part 2:
 DO I=NF2+1,NF
 
   !Part 3:
    ME = IDS(1,I)

	NXX = NX(I)    
	NYY = NY(I)
 
   !Part 4:
    Nu  = WTB(I)/WB(1,I)
    RNu = WTB(I)

   !Part 5:

    DNuXC(ME)  = DNuXC(ME)  + Nu  * NXX
    DNuYC(ME)  = DNuYC(ME)  + Nu  * NYY

    DRNuXC(ME) = DRNuXC(ME) + RNu * NXX
    DRNuYC(ME) = DRNuYC(ME) + RNu * NYY

 End Do
 
!Part 6:
 DO I=NF1+1,NF2
 
   !Part 7:
    ME = IDS(1,I)
	NE = IDS(2,I)

	NXX = NX(I)    
	NYY = NY(I)
 
   !Part 8:
    Nu  = 0.5*( WTNP1(ME)/WNP1(1,ME) + WTNP1(NE)/WNP1(1,NE) )
    RNu = 0.5*( WTNP1(ME) + WTNP1(NE) )

   !Part 9:
    DNuXC(ME)  = DNuXC(ME)  + Nu  * NXX
    DNuYC(ME)  = DNuYC(ME)  + Nu  * NYY

    DNuXC(NE)  = DNuXC(NE)  - Nu  * NXX
    DNuYC(NE)  = DNuYC(NE)  - Nu  * NYY

    DRNuXC(ME) = DRNuXC(ME) + RNu * NXX
    DRNuYC(ME) = DRNuYC(ME) + RNu * NYY

    DRNuXC(NE) = DRNuXC(NE) - RNu * NXX
    DRNuYC(NE) = DRNuYC(NE) - RNu * NYY

 End Do
 
!Part 10: 
 DO I=1,NC

    AA = A(I)

    DNuXC(I)  = DNuXC(I)  / AA
    DNuYC(I)  = DNuYC(I)  / AA

    DRNuXC(I) = DRNuXC(I) / AA
    DRNuYC(I) = DRNuYC(I) / AA
 End Do
!*********************************************************************************************
 End
!###########################################################################################