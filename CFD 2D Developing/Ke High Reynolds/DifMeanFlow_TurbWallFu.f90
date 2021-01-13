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
!// Chief Developer: N. msnkre,Aerospace eng. Amirkabir University of Technology           //!
!// Supervisor: Dr. h. hdhrnuidn,Aerospace eng. Amirkabir University of Technology         //!
!// Date: May.,04,2018                                                                     //!
!// Developed by: N. msnkre,Aerospace Eng.,Amirkabir University of Technology              //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied,Modified and Redistributed for Non-Commercial Use.                    //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DifMeanFlow_TurbWallFu(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,DA,GM,PrL,PrT,NX,NY,MR,Mu,Mut,&
                                WNP1,TurbQ,WB,DUX,DUY,DVX,DVY,DTX,DTY,TauWall,Dif)
 Implicit None
!*******************************************************************************************
INTEGER                           ,INTENT(IN)        ::DIM
INTEGER                           ,INTENT(IN)        ::NC
INTEGER                           ,INTENT(IN)        ::NFW1
INTEGER                           ,INTENT(IN)        ::NFW2
INTEGER                           ,INTENT(IN)        ::NF
INTEGER                           ,INTENT(IN)        ::NF1
INTEGER                           ,INTENT(IN)        ::NF2
INTEGER   ,DIMENSION(1:4,1:DIM)   ,INTENT(IN)        ::IDS
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DA
REAL(8)                           ,INTENT(IN)        ::GM
REAL(8)                           ,INTENT(IN)        ::PRL
REAL(8)                           ,INTENT(IN)        ::PRT
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::NX
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::NY
REAL(8)                           ,INTENT(IN)        ::MR
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::MU
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::MUT
REAL(8)   ,DIMENSION(1:4,1:DIM)   ,INTENT(IN)        ::WNP1
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::TURBQ
REAL(8)   ,DIMENSION(1:5,1:DIM)   ,INTENT(IN)        ::WB
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DUX
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DUY
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DVX
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DVY
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DTX
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::DTY
REAL(8)   ,DIMENSION(1:DIM)       ,INTENT(IN)        ::TAUWALL
REAL(8)   ,DIMENSION(1:4,1:DIM)   ,INTENT(OUT)       ::DIF

INTEGER                                              ::I,ME,NE
REAL(8)                                              ::U,V,NXX,NYY,F1,F2,F3,F4,QX,QY,K,Mum,Txx,Tyy,Txy,MumL,MumT
REAL(8)                                              ::Uii,Rhom,Km,sin,cos,dx,dy,DAA
REAL(8),DIMENSION(1:2,1:Dim)                         ::WTNP1
REAL(8),DIMENSION(1:Dim)                             ::X,Y,XC,YC
!*******************************************************************************************
!Part 1:
 DO I=1,NC
    Dif(2,I) = 0.0
    Dif(3,I) = 0.0
    Dif(4,I) = 0.0
 END DO

!Part 2:
 DO I=NF1+1,NF2
  
   !Part 3:        
    ME = IDS(1,I)        
	NE = IDS(2,I)

   !Part 4:
    NXX = NX(I) 
    NYY = NY(I) 

   !Part 5:
    U    = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )
    V    = 0.5*( WNP1(3,ME)/WNP1(1,ME) + WNP1(3,NE)/WNP1(1,NE) )
    MumL = 0.5*( Mu(ME) + Mu(NE) )
	MumT = 0.5*( Mut(ME) + Mut(NE) )
    Km   = 0.5*( TurbQ(ME) + TurbQ(NE) )
    
   !Part 6:
    Uii = (DUX(I)+DVY(I))/1.5
    Txx = -(MumL+1.0*MumT) * ( 2.0*DUX(I)-Uii    ) + Km/1.5 
	Txy = -(MumL+1.0*MumT) * (   DVX(I)+DUY(I) )   
    Tyy = -(MumL+1.0*MumT) * ( 2.0*DVY(I)-Uii    ) + Km/1.5 

   !Part 7:
    K  = MumL / ((GM-1)*PrL) + MumT / ((GM-1)*PrT)    
    QX =-K*DTX(I)
    QY =-K*DTY(I)
	  
   !Part 8:
    F2 =  Txx*NXX + Txy*NYY
    F3 =  Txy*NXX + Tyy*NYY
    F4 = (U*Txx+V*Txy+QX)*NXX + (U*Txy+V*Tyy+QY)*NYY
    
   !Part 9:
    Dif(2,ME) = Dif(2,ME) + F2
    Dif(3,ME) = Dif(3,ME) + F3
    Dif(4,ME) = Dif(4,ME) + F4
  
   !Part 10:
    Dif(2,NE) = Dif(2,NE) - F2
    Dif(3,NE) = Dif(3,NE) - F3
    Dif(4,NE) = Dif(4,NE) - F4

 END DO

!Part 11:
 DO I=NFW1+1,NFW2
  
   !Part 12:
    ME = IDS(1,I)

   !Part 13:
    NXX = NX(I) 
    NYY = NY(I) 
    DAA = DA(I)

   !Part 14:
	dx  = -NYY
	dy  =  NXX
	sin =  dx/DAA 
	cos =  dy/DAA

   !Part 15:
    Txx = 2*TauWall(I)*sin*cos
	Txy = TauWall(I)*(cos*cos-sin*sin)
    Tyy =-Txx 

   !Part 16:
    F2 =  Txx*NXX + Txy*NYY
    F3 =  Txy*NXX + Tyy*NYY
    
   !Part 17:
    Dif(2,ME) = Dif(2,ME) + F2
    Dif(3,ME) = Dif(3,ME) + F3
        
 End do

!Part 18:
 DO I=NFW2+1,NF
  
   !Part 19:
    ME = IDS(1,I)

   !Part 20:
    NXX = NX(I) 
    NYY = NY(I) 

   !Part 21:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 22:
    Uii = (DUX(I)+DVY(I))/1.5
    MumL = Mu(ME)
	MumT = Mut(ME)

   !Part 23:
    Txx = -(MumL+1.0*MumT) * ( 2.0*DUX(I)-Uii    ) + WTNP1(1,ME)/1.5
	Txy = -(MumL+1.0*MumT) * (   DVX(I)+DUY(I) )
    Tyy = -(MumL+1.0*MumT) * ( 2.0*DVY(I)-Uii    ) + WTNP1(1,ME)/1.5 

   !Part 24:
    K  = MumL / ((GM-1)*PrL) + MumT / ((GM-1)*PrT)    
    QX =-K*DTX(I)
    QY =-K*DTY(I)

   !Part 25:
    F2 =  Txx*NXX + Txy*NYY
    F3 =  Txy*NXX + Tyy*NYY
    F4 = (U*Txx+V*Txy+QX)*NXX + (U*Txy+V*Tyy+QY)*NYY
    
   !Part 26:
    Dif(2,ME) = Dif(2,ME) + F2
    Dif(3,ME) = Dif(3,ME) + F3
    Dif(4,ME) = Dif(4,ME) + F4
        
 End do

!Part 27:
 Do i=1,NC
    Dif(2,I) = MR*Dif(2,I) 
	Dif(3,I) = MR*Dif(3,I) 
	Dif(4,I) = MR*Dif(4,I) 
 End do
!*******************************************************************************************
 End
!###########################################################################################

