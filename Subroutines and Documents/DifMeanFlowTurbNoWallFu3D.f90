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
!// Date: Aug., 30, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DifMeanFlowTurbNoWallFu3D(Dim,NC,NF1,NF2,NFW1,NFW2,NF,IDS,GM,PrL,PrT,MR,Mu,Mut,NX,NY,NZ,DTX,DTY,DTZ,DUX,DUY,DUZ,&
                                      DVX,DVY,DVZ,DWX,DWY,DWZ,WNP1,WB,TurbQ,Dif)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NFW1,NFW2,NF,IDS,GM,PrL,PrT,MR,Mu,Mut,NX,NY,NZ,DTX,DTY,DTZ,DUX,DUY,DUZ,DVX,DVY,&
                DVZ,DWX,DWY,DWZ,WNP1,WB,TurbQ
 Intent(Out  )::Dif

 Integer::Dim,I,j,NC,NF1,NF2,NFW1,NFW2,NF,ME,NE
 Real(8)::U,V,W,NXX,NYY,NZZ,F2,F3,F4,F5,PrL,PrT,QX,QY,QZ,K,GM,MR,MumL,MumT,Mum,Uii,Sxx,Sxy,Sxz,Syx,Syy,&
          Syz,Szx,Szy,Szz,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz,Ro1,Ro2,Ro,Fluc
 Real(8),Dimension(1:Dim)::Mu,NX,NY,NZ,DTX,DTY,DTZ,DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,Mut,TurbQ
 Real(8),Dimension(1:5,1:Dim)::WNP1,Dif
 Real(8),Dimension(1:6,1:Dim)::WB
 Integer,Dimension(1:6,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 DO I=1,NC
    Dif(2,I) = 0.0
    Dif(3,I) = 0.0
    Dif(4,I) = 0.0
    Dif(5,I) = 0.0
 END DO

!Part 2:
 DO I=NF1+1,NF2
  
   !Part 3:        
    ME = IDS(1,I)        
	NE = IDS(2,I)

   !Part 4:
    NXX = NX(I) 
    NYY = NY(I)  
    NZZ = NZ(I)


    !Part 5:
    Ro1 = WNP1(1,ME)
	Ro2 = WNP1(1,NE)

    U  = 0.5 * ( WNP1(2,ME)/Ro1 + WNP1(2,NE)/Ro2 )
    V  = 0.5 * ( WNP1(3,ME)/Ro1 + WNP1(3,NE)/Ro2 )
    W  = 0.5 * ( WNP1(4,ME)/Ro1 + WNP1(4,NE)/Ro2 )

   
    MumL = 0.5*(Mu(ME) +Mu(NE) )
	MumT = 0.5*(Mut(ME)+Mut(NE))
    Mum = MumL+MumT

	Fluc = 0.5*( TurbQ(ME)+TurbQ(NE) ) / 1.5
    
   !Part 6:
    Uii = (DUX(I)+DVY(I)+DWZ(I))/1.5

    Sxx = 2*DUX(I)-Uii ; Sxy = DUY(I)+DVX(I) ; Sxz = DUZ(I)+DWX(I)
                         Syy = 2*DVY(I)-Uii  ; Syz = DVZ(I)+DWY(I)
                                               Szz = 2*DWZ(I)-Uii
    
    
    Txx = -Mum * Sxx + Fluc ; Txy = -Mum * Sxy        ; Txz = -Mum * Sxz
    Tyx = Txy               ; Tyy = -Mum * Syy + Fluc ; Tyz = -Mum * Syz
    Tzx = Txz               ; Tzy = Tyz               ; Tzz = -Mum * Szz + Fluc
    
   !Part 7:	
    K  = MumL / ((GM-1)*PrL) + MumT / ((GM-1)*PrT)     
    Qx =-K*DTX(I)
    Qy =-K*DTY(I)
    Qz =-K*DTZ(I)

   !Part 8:
    F2 =  Txx*NXX + Txy*NYY + Txz*Nzz
    F3 =  Tyx*NXX + Tyy*NYY + Tyz*Nzz
    F4 =  Tzx*NXX + Tzy*NYY + Tzz*Nzz
    F5 = (U*Txx+V*Tyx+W*Tzx+Qx)*Nxx + (U*Txy+V*Tyy+W*Tzy+Qy)*Nyy + (U*Txz+V*Tyz+W*Tzz+Qz)*Nzz
         
   !Part 9:
    Dif(2,ME) = Dif(2,ME) + F2
    Dif(3,ME) = Dif(3,ME) + F3
    Dif(4,ME) = Dif(4,ME) + F4
    Dif(5,ME) = Dif(5,ME) + F5
  
   !Part 10:
    Dif(2,NE) = Dif(2,NE) - F2
    Dif(3,NE) = Dif(3,NE) - F3
    Dif(4,NE) = Dif(4,NE) - F4
    Dif(5,NE) = Dif(5,NE) - F5


 END DO

!Part 11:
 DO I=NFW1+1,NFW2
  
   !Part 12:
    ME = IDS(1,I) 

   !Part 13:
    NXX = NX(I) 
    NYY = NY(I)  
    NZZ = NZ(I)

   !Part 14:
    Mum = Mu(ME)

   !Part 15:
    Uii = (DUX(I)+DVY(I)+DWZ(I))/1.5

    Sxx = 2*DUX(I)-Uii ; Sxy = DUY(I)+DVX(I) ; Sxz = DUZ(I)+DWX(I)
                         Syy = 2*DVY(I)-Uii  ; Syz = DVZ(I)+DWY(I)
                                               Szz = 2*DWZ(I)-Uii

    Txx = -Mum * Sxx ; Txy = -Mum * Sxy  ; Txz = -Mum * Sxz
    Tyx = Txy        ; Tyy = -Mum * Syy  ; Tyz = -Mum * Syz
    Tzx = Txz        ; Tzy = Tyz         ; Tzz = -Mum * Szz 
   	  
   !Part 16:
    F2 =  Txx*NXX + Txy*NYY + Txz*Nzz
    F3 =  Tyx*NXX + Tyy*NYY + Tyz*Nzz
    F4 =  Tzx*NXX + Tzy*NYY + Tzz*Nzz
    
    Dif(2,ME) = Dif(2,ME) + F2
    Dif(3,ME) = Dif(3,ME) + F3
    Dif(4,ME) = Dif(4,ME) + F4

        
 END DO

 
!Part 17:
 DO I=NFW2+1,NF
  
   !Part 18:
    ME = IDS(1,I) 

   !Part 19:
    NXX = NX(I) 
    NYY = NY(I)  
    NZZ = NZ(I)

   !Part 20:
    MumL = Mu(ME)
    MumT = Mut(ME)
    Mum = MumL+MumT

    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)

   !Part 21:
    Uii = (DUX(I)+DVY(I)+DWZ(I))/1.5

    Sxx = 2*DUX(I)-Uii ; Sxy = DUY(I)+DVX(I) ; Sxz = DUZ(I)+DWX(I)
                         Syy = 2*DVY(I)-Uii  ; Syz = DVZ(I)+DWY(I)
                                               Szz = 2*DWZ(I)-Uii
    
    Txx = -Mum * Sxx ; Txy = -Mum * Sxy ; Txz = -Mum * Sxz
    Tyx = Txy        ; Tyy = -Mum * Syy ; Tyz = -Mum * Syz
    Tzx = Txz        ; Tzy = Tyz        ; Tzz = -Mum * Szz
    
   !Part 22:	
    K  = MumL / ((GM-1)*PrL) + MumT / ((GM-1)*PrT)     
    Qx =-K*DTX(I)
    Qy =-K*DTY(I)
    Qz =-K*DTZ(I)

   !Part 23:
    F2 =  Txx*NXX + Txy*NYY + Txz*Nzz
    F3 =  Tyx*NXX + Tyy*NYY + Tyz*Nzz
    F4 =  Tzx*NXX + Tzy*NYY + Tzz*Nzz
    F5 = (U*Txx+V*Tyx+W*Tzx+Qx)*Nxx + (U*Txy+V*Tyy+W*Tzy+Qy)*Nyy + (U*Txz+V*Tyz+W*Tzz+Qz)*Nzz
    
   !Part 24:
    Dif(2,ME) = Dif(2,ME) + F2
    Dif(3,ME) = Dif(3,ME) + F3
    Dif(4,ME) = Dif(4,ME) + F4
    Dif(5,ME) = Dif(5,ME) + F5

        
 END DO

!Part 25:
 Do i=1,NC
    Dif(2,I) = MR*Dif(2,I) 
	Dif(3,I) = MR*Dif(3,I) 
	Dif(4,I) = MR*Dif(4,I) 
	Dif(5,I) = MR*Dif(5,I) 
 End do
!*********************************************************************************************
 End
!###########################################################################################




!!!!Part 12:
!!! DO I=NFW2+1,0  !NFF
!!!  
!!!   !Part 13:
!!!    ME = IDS(1,I) 
!!!
!!!   !Part 4:
!!!    NXX = NX(I) 
!!!    NYY = NY(I)  
!!!    NZZ = NZ(I)
!!!
!!!   !Part 5:
!!!    Ro = WB(1,I)
!!!    U  = WB(2,I)/Ro
!!!    V  = WB(3,I)/Ro
!!!    W  = WB(4,I)/Ro
!!!    MumL = Mu(ME)
!!!
!!!   !Part 6:
!!!    Uii = (DUX(I)+DVY(I)+DWZ(I))/1.5
!!!
!!!    Sxx =2*DUX(I)-Uii      ;   Sxy =  DUY(I)+DVX(I)   ;   Sxz =  DUZ(I)+DWX(I)
!!!    Syx =  DVX(I)+DUY(I)   ;   Syy =2*DVY(I)-Uii      ;   Syz =  DVZ(I)+DWY(I)
!!!    Szx =  DWX(I)+DUZ(I)   ;   Szy =  DWY(I)+DVZ(I)   ;   Szz =2*DWZ(I)-Uii
!!!
!!!    Txx = -MumL * Sxx  ;   Txy = -MumL * Sxy  ;   Txz = -MumL * Sxz
!!!    Tyx = -MumL * Syx  ;   Tyy = -MumL * Syy  ;   Tyz = -MumL * Syz
!!!    Tzx = -MumL * Szx  ;   Tzy = -MumL * Szy  ;   Tzz = -MumL * Szz
!!!   
!!!    K  = MumL / ((GM-1)*PrL)    
!!!    Qx =-K*DTX(I)
!!!    Qy =-K*DTY(I)
!!!    Qz =-K*DTZ(I)
!!!
!!!   !Part 8:
!!!    F2 =  Txx*NXX + Txy*NYY + Txz*Nzz
!!!    F3 =  Tyx*NXX + Tyy*NYY + Tyz*Nzz
!!!    F4 =  Tzx*NXX + Tzy*NYY + Tzz*Nzz
!!!    F5 = (U*Txx+V*Tyx+W*Tzx+Qx)*Nxx + (U*Txy+V*Tyy+W*Tzy+Qy)*Nyy + (U*Txz+V*Tyz+W*Tzz+Qz)*Nzz
!!!
!!!   !Part 9:
!!!    Dif(2,ME) = Dif(2,ME) + F2
!!!    Dif(3,ME) = Dif(3,ME) + F3
!!!    Dif(4,ME) = Dif(4,ME) + F4
!!!
!!!        
!!! END DO

