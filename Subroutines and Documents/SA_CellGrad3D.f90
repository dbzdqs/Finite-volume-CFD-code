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
 Subroutine SA_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WTNP1,WTB,WB,DNuX_C,DNuY_C,DNuZ_C,DRNuX_C,DRNuY_C,DRNuZ_C)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WTNP1,WTB,WB
 Intent(Out  )::DNuX_C,DNuY_C,DNuZ_C,DRNuX_C,DRNuY_C,DRNuZ_C

 Integer::Dim,I,NC,NF,ME,NE,NF1,NF2
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8)::Nu,RNu,Rho,Volume,NXX,NYY,NZZ
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,Vol,DNuX_C,DNuY_C,DNuZ_C,DRNuX_C,DRNuY_C,DRNuZ_C
!********************************************************************************************* 
!Part 1:
 DO I=1,NC
    DNuX_C(I) = 0.0
    DNuY_C(I) = 0.0
    DNuZ_C(I) = 0.0
    
    DRNuX_C(I) = 0.0
    DRNuY_C(I) = 0.0
    DRNuZ_C(I) = 0.0
 End Do

!Part 2:
 DO I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 4:
    Nu  = 0.5*( WTNP1(1,ME)/WNP1(1,ME)+WTNP1(1,NE)/WNP1(1,NE) )
    RNu = 0.5*( WTNP1(2,ME)+WTNP1(2,NE) )

   !Part 5:
    NXX = NX(I)
    NYY = NY(I)
    NZZ = NZ(I)

    DNuX_C(ME)  = DNuX_C(ME)  + Nu*NXX
    DNuY_C(ME)  = DNuY_C(ME)  + Nu*NYY
    DNuZ_C(ME)  = DNuZ_C(ME)  + Nu*NZZ
    
    DRNuX_C(ME) = DRNuX_C(ME) + RNu*NXX
    DRNuY_C(ME) = DRNuY_C(ME) + RNu*NYY
    DRNuZ_C(ME) = DRNuZ_C(ME) + RNu*NZZ

    DNuX_C(NE)  = DNuX_C(NE)  - Nu*NXX
    DNuY_C(NE)  = DNuY_C(NE)  - Nu*NYY
    DNuZ_C(NE)  = DNuZ_C(NE)  - Nu*NZZ
    
    DRNuX_C(NE) = DRNuX_C(NE) - RNu*NXX
    DRNuY_C(NE) = DRNuY_C(NE) - RNu*NYY
    DRNuZ_C(NE) = DRNuZ_C(NE) - RNu*NZZ
 
 End Do

!Part 6:
 DO I=NF2+1,NF

   !Part 7:
    ME = IDS(1,I)

    NXX = NX(I)
    NYY = NY(I)
    NZZ = NZ(I)

   !Part 8:
    Rho = WB(1,I)
    Nu  = WTB(1,I)/Rho
    RNu = WTB(2,I)/Rho

   !Part 9:
    DNuX_C(ME)  = DNuX_C(ME)  + Nu*NXX
    DNuY_C(ME)  = DNuY_C(ME)  + Nu*NYY
    DNuZ_C(ME)  = DNuZ_C(ME)  + Nu*NZZ
    
    DRNuX_C(ME) = DRNuX_C(ME) + RNu*NXX
    DRNuY_C(ME) = DRNuY_C(ME) + RNu*NYY
    DRNuZ_C(ME) = DRNuZ_C(ME) + RNu*NZZ

 End Do 
 
!Part 10:
 DO I=1,NC
    Volume = Vol(I)

    DNuX_C(I)  = DNuX_C(I)  / Volume
    DNuY_C(I)  = DNuY_C(I)  / Volume
    DNuZ_C(I)  = DNuZ_C(I)  / Volume

    DRNuX_C(I) = DRNuX_C(I) / Volume
    DRNuY_C(I) = DRNuY_C(I) / Volume
    DRNuZ_C(I) = DRNuZ_C(I) / Volume
 End Do
!*********************************************************************************************
 End
!###########################################################################################

