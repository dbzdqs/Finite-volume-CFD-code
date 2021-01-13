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
 Subroutine Kw_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WTNP1,WTB,WB,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WTNP1,WTB,WB
 Intent(Out  )::DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C

 Integer::Dim,I,NC,NF,ME,NE,NF1,NF2
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8)::K,Omeg,Rho,Volume,NXX,NYY,NZZ
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,Vol,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C
!********************************************************************************************* 
!Part 1:
 DO I=1,NC
    DKX_C(I) = 0.0
    DKY_C(I) = 0.0
    DKZ_C(I) = 0.0
    
    DOmegX_C(I) = 0.0
    DOmegY_C(I) = 0.0
    DOmegZ_C(I) = 0.0
 End Do

!Part 2:
 DO I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 4:
    K    = 0.5*( WTNP1(1,ME)/WNP1(1,ME)+WTNP1(1,NE)/WNP1(1,NE) )
    Omeg = 0.5*( WTNP1(2,ME)/WNP1(1,ME)+WTNP1(2,NE)/WNP1(1,NE) )

   !Part 5:
    NXX = NX(I)
    NYY = NY(I)
    NZZ = NZ(I)

    DKX_C(ME)    = DKX_C(ME)    + K*NXX
    DKY_C(ME)    = DKY_C(ME)    + K*NYY
    DKZ_C(ME)    = DKZ_C(ME)    + K*NZZ
    
    DOmegX_C(ME) = DOmegX_C(ME) + Omeg*NXX
    DOmegY_C(ME) = DOmegY_C(ME) + Omeg*NYY
    DOmegZ_C(ME) = DOmegZ_C(ME) + Omeg*NZZ

    DKX_C(NE)    = DKX_C(NE)    - K*NXX
    DKY_C(NE)    = DKY_C(NE)    - K*NYY
    DKZ_C(NE)    = DKZ_C(NE)    - K*NZZ
    
    DOmegX_C(NE) = DOmegX_C(NE) - Omeg*NXX
    DOmegY_C(NE) = DOmegY_C(NE) - Omeg*NYY
    DOmegZ_C(NE) = DOmegZ_C(NE) - Omeg*NZZ
 
 End Do

!Part 6:
 DO I=NF2+1,NF

   !Part 7:
    ME = IDS(1,I)

    NXX = NX(I)
    NYY = NY(I)
    NZZ = NZ(I)

   !Part 8:
    Rho  = WB(1,I)
    K    = WTB(1,I)/Rho
    Omeg = WTB(2,I)/Rho

   !Part 9:
    DKX_C(ME)    = DKX_C(ME)    + K*NXX
    DKY_C(ME)    = DKY_C(ME)    + K*NYY
    DKZ_C(ME)    = DKZ_C(ME)    + K*NZZ
    
    DOmegX_C(ME) = DOmegX_C(ME) + Omeg*NXX
    DOmegY_C(ME) = DOmegY_C(ME) + Omeg*NYY
    DOmegZ_C(ME) = DOmegZ_C(ME) + Omeg*NZZ

 End Do 
 
!Part 10:
 DO I=1,NC
    Volume = Vol(I)

    DKX_C(I)    = DKX_C(I)   / Volume
    DKY_C(I)    = DKY_C(I)   / Volume
    DKZ_C(I)    = DKZ_C(I)   / Volume

    DOmegX_C(I) = DOmegX_C(I) / Volume
    DOmegY_C(I) = DOmegY_C(I) / Volume
    DOmegZ_C(I) = DOmegZ_C(I) / Volume
 End Do
!*********************************************************************************************
 End
!###########################################################################################

