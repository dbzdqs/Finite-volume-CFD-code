!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the gradient at Cell Center of 2 eq. Turbulence Model 2D      //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F055F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Velocity_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB
 Intent(Out  )::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C

 Integer::Dim,I,NC,NF,ME,NE,NF1,NF2
 Real(8)::U,V,W,Rho,Volume,NXX,NYY,NZZ
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,NZ,Vol,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
!******************************************************************************************* 
!Part 1:
 DO I=1,NC
    DUX_C(I) = 0.0
    DUY_C(I) = 0.0
    DUZ_C(I) = 0.0
    
    DVX_C(I) = 0.0
    DVY_C(I) = 0.0
    DVZ_C(I) = 0.0
    
    DWX_C(I) = 0.0
    DWY_C(I) = 0.0
    DWZ_C(I) = 0.0
 End Do

!Part 2:
 DO I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 4:    
    U  = 0.5*( WNP1(2,ME)/WNP1(1,ME)+WNP1(2,NE)/WNP1(1,NE) )
    V  = 0.5*( WNP1(3,ME)/WNP1(1,ME)+WNP1(3,NE)/WNP1(1,NE) )
    W  = 0.5*( WNP1(4,ME)/WNP1(1,ME)+WNP1(4,NE)/WNP1(1,NE) )

   !Part 5:
    NXX = NX(I)
    NYY = NY(I)
    NZZ = NZ(I)

    DUX_C(ME)  = DUX_C(ME)  + U*NXX
    DUY_C(ME)  = DUY_C(ME)  + U*NYY
    DUZ_C(ME)  = DUZ_C(ME)  + U*NZZ
    
    DVX_C(ME)  = DVX_C(ME)  + V*NXX
    DVY_C(ME)  = DVY_C(ME)  + V*NYY
    DVZ_C(ME)  = DVZ_C(ME)  + V*NZZ
    
    DWX_C(ME)  = DWX_C(ME)  + W*NXX
    DWY_C(ME)  = DWY_C(ME)  + W*NYY
    DWZ_C(ME)  = DWZ_C(ME)  + W*NZZ
    
    
    DUX_C(NE)  = DUX_C(NE)  - U*NXX
    DUY_C(NE)  = DUY_C(NE)  - U*NYY
    DUZ_C(NE)  = DUZ_C(NE)  - U*NZZ
    
    DVX_C(NE)  = DVX_C(NE)  - V*NXX
    DVY_C(NE)  = DVY_C(NE)  - V*NYY
    DVZ_C(NE)  = DVZ_C(NE)  - V*NZZ
    
    DWX_C(NE)  = DWX_C(NE)  - W*NXX
    DWY_C(NE)  = DWY_C(NE)  - W*NYY
    DWZ_C(NE)  = DWZ_C(NE)  - W*NZZ

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
    U   = WB(2,I) /Rho
    V   = WB(3,I) /Rho
    W   = WB(4,I) /Rho

   !Part 9:
    DUX_C(ME)  = DUX_C(ME)  + U*NXX
    DUY_C(ME)  = DUY_C(ME)  + U*NYY
    DUZ_C(ME)  = DUZ_C(ME)  + U*NZZ
    
    DVX_C(ME)  = DVX_C(ME)  + V*NXX
    DVY_C(ME)  = DVY_C(ME)  + V*NYY
    DVZ_C(ME)  = DVZ_C(ME)  + V*NZZ
    
    DWX_C(ME)  = DWX_C(ME)  + W*NXX
    DWY_C(ME)  = DWY_C(ME)  + W*NYY
    DWZ_C(ME)  = DWZ_C(ME)  + W*NZZ
 End Do 
 
!Part 10:
 DO I=1,NC
    Volume = Vol(I)

    DUX_C(I)  = DUX_C(I) / Volume 
    DUY_C(I)  = DUY_C(I) / Volume 
    DUZ_C(I)  = DUZ_C(I) / Volume 
    
    DVX_C(I)  = DVX_C(I) / Volume 
    DVY_C(I)  = DVY_C(I) / Volume 
    DVZ_C(I)  = DVZ_C(I) / Volume 
    
    DWX_C(I)  = DWX_C(I) / Volume 
    DWY_C(I)  = DWY_C(I) / Volume 
    DWZ_C(I)  = DWZ_C(I) / Volume 
 End Do 
!*******************************************************************************************
 End
!###########################################################################################

