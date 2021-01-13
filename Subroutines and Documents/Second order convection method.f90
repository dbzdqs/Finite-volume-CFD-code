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
 Subroutine ConMeanFlow_AUSM_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con


 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,U_L,V_L,P_L,E_L,RHO_R,U_R,V_R,P_R,E_R
 Real(8)::PR,PHI,DX,DY

!*********************************************************************************************	

!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA
    
    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       
       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)

      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
    

   !Part 11:
    a_L  = Dsqrt(GM*P_L/Rho_L)
	a_R  = Dsqrt(GM*P_R/Rho_R)

	M_L = (Rho_L*U_L*NXX+Rho_L*V_L*NYY) / (Rho_L * a_L)
	M_R = (Rho_R*U_R*NXX+Rho_R*V_R*NYY) / (Rho_R * a_R)

   !Part 12:
	IF(Dabs(M_L)>1.)Then
	 M_Plus = 0.5*(M_L+Dabs(M_L))
	 P_Plus = 0.5*(M_L+Dabs(M_L))/(M_L)
	Else
     M_Plus = 0.25*(M_L+1.)*(M_L+1.)
	 P_Plus = 0.25*(M_L+1.)*(M_L+1.)*(2.-M_L)
	End If

   !Part 13:
	IF(Dabs(M_R)>1.)Then
	 M_Minus = 0.5*(M_R-Dabs(M_R))
	 P_Minus = 0.5*(M_R-Dabs(M_R))/(M_R)
	Else
     M_Minus =-0.25*(M_R-1.)*(M_R-1.)
	 P_Minus = 0.25*(M_R-1.)*(M_R-1.)*(2.+M_R)
	End If

   !Part 14:
    Mm = M_Plus+M_Minus
	Pm = P_L*P_Plus + P_R*P_Minus 

   !Part 15:
	If(Mm<=0.)Then
     Ro = Rho_R       * a_R
	 RU = Rho_R*U_R       * a_R
	 RV = Rho_R*V_R       * a_R
	 RH =(Rho_R*E_R+P_R) * a_R
	Else
     Ro = Rho_L       * a_L
	 RU = Rho_L*U_L       * a_L
	 RV = Rho_L*V_L       * a_L
	 RH =(Rho_L*E_L+P_L) * a_L
	Endif

   !Part 16:
    F1 = ( Mm * Ro          ) * DAA
    F2 = ( Mm * RU + Pm*NXX ) * DAA
    F3 = ( Mm * RV + Pm*NYY ) * DAA
    F4 = ( Mm * RH          ) * DAA

   !Part 17:
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4

   !Part 18:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    
    
    
    
 End Do
!*********************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Convection Terms of 2D Mean Flow Equations Using AUSM+    //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: Ali Reza Rezaei                                                        //!
!// Developed by: N. msnkre, A. Rezayi Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC2F078F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_AUSM_Plus_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,Con,X,Y)
 Implicit none
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,R,L
 Real(8)::Alpha,Beta,Const1,GM,U,V,Q,Pm,Rho_L,Rho_R,U_L,U_R,V_L,V_R,P_L,P_R,a_L,a_R,astar_L,&
          astar_R,Utot_L,Utot_R,Ucontra_L,Ucontra_R,a_Face,Ma_L,Ma_R,Ma_Plus_L,Ma_minus_R,Phi,&
          Ma_face,MassFlux,P_plus_L,P_minus_R,P_face,H_L,H_R,W,F1,F2,F3,F4,DAA,NXX,NYY,Ro,Pr,DX,DY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,DA,P,Xc,Yc
 Real(8),Dimension(1:4,1:Dim)::Con,WNP1,Limit
 Real(8),Dimension(1:5,1:Dim)::WB    
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 
 
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
!*********************************************************************************************
!Part1 :
 Alpha = ( 3.0 / 16.0 )
 Beta = 0.125
 Const1 = 2.0 * (GM - 1.0 ) / ( GM + 1.0)


!Part 2:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 3:
 DO I=NF2+1,NF

   !Part 4:
    ME = IDS(1,I)

   !Part 5:
    U = WB(2,I) / WB(1,I)
    V = WB(3,I) / WB(1,I)

   !Part 6:
    Q = U*NX(I) + V*NY(I)       
    Pm = WB(5,I)

   !Part 7:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I) 
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 8:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4

 End Do

!Part 9:
 DO I=NF1+1,NF2

   !Part 10:
    L = IDS(1,I)
    R = IDS(2,I)

    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr

      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr

      !Part 13:
       Utot_L  = U_L*U_L + V_L*V_L
       H_L     = ( (  GM * P_L / Rho_L  ) / ( GM - 1.0 ) ) + ( 0.5 * Utot_L )
       astar_L = Dsqrt(H_L  * Const1)

       Utot_R  =  U_R*U_R + V_R*V_R
       H_R     = ( ( GM * P_R / Rho_R ) / ( GM - 1.0 ) ) + ( 0.5 * Utot_R ) 
       astar_R = Dsqrt(H_R  * Const1)   
   
      !Part 14:
       Ucontra_L = U_L*NXX + V_L*NYY
       Ucontra_R = U_R*NXX + V_R*NYY

      !Part 15:
       a_L = ( astar_L*astar_L ) / DMax1( astar_L , Ucontra_L )
       a_R = ( astar_R*astar_R ) / DMax1( astar_R ,-Ucontra_R )

       a_face = DMax1( a_L , a_R )

      !Part 16:    
       Ma_L = Ucontra_L / a_face
       Ma_R = Ucontra_R / a_face    
      
      !Part 17: 
       IF ( ABS(Ma_L) > 1.0 ) Then
        Ma_plus_L = 0.5 * ( Ma_L + ABS(Ma_L) )
       Else 
        Ma_plus_L =  0.25 * ( ( Ma_L + 1.0 )**2.0 ) + Beta * ( ( (Ma_L**2.0) -1 ) ** 2.0 ) 
       End IF

       IF ( ABS(Ma_R) > 1.0 ) Then
        Ma_minus_R = 0.5 * ( Ma_R - ABS(Ma_R) )
       Else 
        Ma_minus_R = -0.25 * ( ( Ma_R - 1.0 )**2.0 )  - Beta * ( ( (Ma_R**2.0) -1 ) ** 2.0 )  
       End IF

      !Part 18:  
       Ma_face = Ma_plus_L + Ma_minus_R       
  
      !Part 19: 
       IF ( ABS(Ma_L) > 1.0 ) Then
        P_plus_L = 0.5* ( 1 + Ma_L/ABS( Ma_L )  )
       Else 
        P_plus_L = 0.25 * ( ( Ma_L+1.0 ) ** 2.0 ) * ( 2.0 - Ma_L )&
                + Alpha * Ma_L* ( ( Ma_L** 2.0 - 1.0 ) ** 2.0 )
       End IF
   
       IF ( ABS(Ma_R) > 1.0 )Then
        P_minus_R = 0.5* ( 1 - Ma_R/ABS( Ma_R )  )
       Else 
        P_minus_R =  0.25 * ( ( Ma_R-1.0 ) ** 2.0 ) * ( 2.0 + Ma_R )&
                - Alpha * Ma_R * ( ( Ma_R** 2.0 - 1.0 ) ** 2.0 )  
       End IF        

      !Part 20:
       P_Face = P_plus_L * P_L + P_minus_R * P_R

      !Part 21:    
       IF ( Ma_face >= 0.0 )Then 
        MassFlux = a_face * Ma_face * Rho_L
       Else
        MassFlux = a_face * Ma_face * Rho_R
       End IF  
    
      !Part 22:
       IF ( MassFlux >= 0.0 )Then

 
        F1 = ( MassFlux                      ) * DAA
        F2 = ( MassFlux * U_L + P_Face * NXX ) * DAA
        F3 = ( MassFlux * V_L + P_Face * NYY ) * DAA
        F4 = ( MassFlux * H_L                ) * DAA
       Else
        F1 = ( MassFlux                      ) * DAA
        F2 = ( MassFlux * U_R + P_Face * NXX ) * DAA
        F3 = ( MassFlux * V_R + P_Face * NYY ) * DAA
        F4 = ( MassFlux * H_R                ) * DAA
       End IF   
   
      !Part 23:
       Con(1,L) = Con(1,L) + F1
       Con(2,L) = Con(2,L) + F2
       Con(3,L) = Con(3,L) + F3
       Con(4,L) = Con(4,L) + F4 

       Con(1,R) = Con(1,R) - F1 
       Con(2,R) = Con(2,R) - F2 
       Con(3,R) = Con(3,R) - F3 
       Con(4,R) = Con(4,R) - F4     

 End Do
!*********************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:  Calculate the Convection Terms of 2D Mean Flow Equations Using AUSM+_Up//!
!//                                                                                      //!
!// Version:                                                                             //!
!// Date:                                                                                //!
!// Developed by: Ali Reza Rezaei                                                        //!
!// Developed by: N. msnkre, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F064F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_AUSM_PlusUP_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P,GWNP1,Xc,Yc,Limit,Con,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::Alpha,Beta,Ku,Kp,Const1,GM,Sigma,Temp1,Temp2,Temp3,U,V,Q,Pm,Rho_L,Rho_R,u_L,u_R,&
          v_L,v_R,P_L,P_R,a_L,a_R,Minf,astar_L,astar_R,Utot_L,Utot_R,Ucontra_L,Ucontra_R,&
          a_Face,Ma_L,Ma_R,Ma_Plus_L,Ma_minus_R,Ma_bar,Mo,M_Co,fa,Ma_p,Ma_face,MassFlux,&
          P_plus_L,P_minus_R,P_u,P_face,H_L,H_R,F1,F2,F3,F4,NXX,NYY,DAA,Ro,Pr,DX,DY,Phi
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,DA,P,Xc,Yc
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con,Limit
 Real(8),Dimension(1:5,1:Dim)::WB    
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:2)::wegt
 Real(8),Dimension(1:2,1:Dim)::XF,YF   
 
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
!*********************************************************************************************	
!Part1 :
 Beta   = 0.125
 Ku     = 0.75
 Kp     = 0.25
 Sigma  = 1.0
 Const1 = 2.0 * (GM - 1.0 ) / ( GM + 1.0)
 

!Part 2:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do  
   
!Part 3:
 DO I=NF2+1,NF

   !Part 4:
    ME = IDS(1,I)

   !Part 5:
    U = WB(2,I) / WB(1,I)
    V = WB(3,I) / WB(1,I)

   !Part 6:
    Q  = U*NX(I) + V*NY(I)       
    Pm = WB(5,I)

   !Part 7:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I) 
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 8:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4

 End Do

!Part 9:
 DO I=NF1+1,NF2

   !Part 10:
    L = IDS(1,I)
    R = IDS(2,I)

    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr

      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr

      !Part 13:
      !Left
       Utot_L  = u_L*u_L + v_L*v_L
       H_L     = ( (GM * P_L / Rho_L) / (GM - 1.0) ) + (0.5 * Utot_L)
       astar_L = Dsqrt(H_L  * Const1)
     
      !Right
       Utot_R  = u_R*u_R + v_R*v_R    
       H_R     = ( (GM * P_R / Rho_R) / (GM - 1.0) ) + (0.5 * Utot_R)
       astar_R = Dsqrt(H_R  * Const1)    

      !Part 14:
       Ucontra_L = u_L*NXX + v_L*NYY  
       Ucontra_R = u_R*NXX + v_R*NYY
    
      !Part 15:
      !Left
       a_L = ( astar_L*astar_L ) / DMax1( astar_L , Ucontra_L )
       a_R = ( astar_R*astar_R ) / DMax1( astar_R ,-Ucontra_R )

       a_face = DMax1( a_L , a_R )

      !Part 16:    
       Ma_L = Ucontra_L / a_face
       Ma_R = Ucontra_R / a_face    

       !Part 17:    
       IF( ABS(Ma_L) >= 1.0 )Then
        Ma_plus_L = 0.5 * ( Ma_L + ABS(Ma_L) )
       Else 
        Ma_plus_L =  0.25 * ( ( Ma_L + 1.0 )**2.0 ) + Beta * ( ( (Ma_L**2.0) -1 ) ** 2.0 ) 
       End IF

       IF( ABS(Ma_R) >= 1.0 )Then
        Ma_minus_R = 0.5 * ( Ma_R - ABS(Ma_R) )
       Else 
        Ma_minus_R = -0.25 * ( ( Ma_R - 1.0 )**2.0 )  - Beta * ( ( (Ma_R**2.0) -1 ) ** 2.0 )  
       End IF    
       
      !Part 18:
       Ma_bar = 0.5 * ( Ma_L*Ma_L + Ma_R*Ma_R )   
   
      !Part 19:
       M_Co = Minf*Minf
       Mo = min(1.0,max(Ma_bar,M_Co))
       Mo = Dsqrt(Mo)   
       fa = Mo * ( 2.0 - Mo )

      !Part 20: 
       Temp1 = 1.0 - (Sigma * Ma_bar )
       Temp2 = 0.5 * ( Rho_L + Rho_R ) * ( a_face *a_face )
       Ma_p = -(Kp/fa) * Max(Temp1,0.0) * ( P_R - P_L ) / Temp2        

      !Part 21:
       Ma_face = Ma_plus_L + Ma_minus_R + Ma_p     

      !Part 22:    
       IF( Ma_face > 0.0 )Then 
        MassFlux = a_face * Ma_face * Rho_L
       Else
        MassFlux = a_face * Ma_face * Rho_R
       End IF         

      !Part 23:
       Alpha = ( 3.0 / 16.0 ) * ( -4.0 + 5.0 * fa * fa ) 
 
       IF( ABS(Ma_L) >= 1.0 )Then
        P_plus_L = 0.5* ( 1 + Ma_L/ABS( Ma_L )  )
       Else 
        P_plus_L = 0.25 * ( ( Ma_L+1.0 ) ** 2.0 ) * ( 2.0 - Ma_L ) + &
                Alpha * Ma_L* ( ( Ma_L** 2.0 - 1.0 ) ** 2.0 )
       End IF
   
       IF( ABS(Ma_R) >= 1.0 )Then
        P_minus_R = 0.5* ( 1 - Ma_R/ABS( Ma_R )  )
       Else 
        P_minus_R = 0.25 * ( ( Ma_R-1.0 ) ** 2.0 ) * ( 2.0 + Ma_R ) - &
                    Alpha * Ma_R * ( ( Ma_R** 2.0 - 1.0 ) ** 2.0 )  
       End IF        

      !Part 24:
       P_u = - Ku*P_plus_L*P_minus_R*( Rho_L + Rho_R )*fa *a_Face*( Ucontra_R - Ucontra_L )    
      
      !Part 25:
       P_Face = P_plus_L * P_L + P_minus_R * P_R + P_u    

      !Part 26:
       IF( Ma_face >= 0. )Then 
        F1 = ( MassFlux                      ) * DAA
        F2 = ( MassFlux * u_L + P_Face * NXX ) * DAA
        F3 = ( MassFlux * v_L + P_Face * NYY ) * DAA
        F4 = ( MassFlux * H_L                ) * DAA
       Else
        F1 = ( MassFlux                      ) * DAA 
        F2 = ( MassFlux * u_R + P_Face * NXX ) * DAA
        F3 = ( MassFlux * v_R + P_Face * NYY ) * DAA
        F4 = ( MassFlux * H_R                ) * DAA
       End IF          

      !Part 27:
       Con(1,L) = Con(1,L) + F1
       Con(2,L) = Con(2,L) + F2
       Con(3,L) = Con(3,L) + F3
       Con(4,L) = Con(4,L) + F4
	    
       Con(1,R) = Con(1,R) - F1
       Con(2,R) = Con(2,R) - F2
       Con(3,R) = Con(3,R) - F3
       Con(4,R) = Con(4,R) - F4   

    
 End Do 
!********************************************************************************************* 
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations Using ROE METHOD//!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: MARCH, 05, 2016                                                                //!
!// Developed by: N. msnkre, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F002F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ZHA_CUSP2_METHOD_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L,P_L,E_L
 Real(8)::R_R,U_R,V_R,H_R,P_R,E_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,DEL_UC,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 Real(8):: R_BAR,A_BAR,P_BAR,U_BAR,V_BAR,E_BAR,H_BAR,UC_BAR
 Real(8)::AC_L,AC_R,AC_BAR,M_BAR
 Real(8)::DLT_P,DLT_M,M_LP,M_LM,M_RP,M_RM
 Real(8)::BET_L,BET_R,M_HALF,M_HALFP,M_HALFM
 Real(8)::ALF_PL,ALF_MR,AC_P,AC_M,P_PL,P_MR,S_PL,S_MR,D_PL,D_MR
 Real(8)::F1_C,F2_C,F3_C,F4_C
 Real(8)::F1_P,F2_P,F3_P,F4_P
 
 Real(8)::MM_L,MM_R
 Real(8)::ALF_L,ALF_R,UP_L,UM_R,RU_HALF,ALF
 INTEGER::IZHA
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,RHO_R
 Real(8)::PR,PHI,DX,DY
 
 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)
    
    
   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)
    
    !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA
    
    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)
       R_L = Rho_L
       H_L = (R_L*E_L+P_L)/R_L
       A_L  = DSQRT(GM*P_L/R_L)

       
      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
       R_R = Rho_R
       H_R = (R_R*E_R+P_R)/R_R
       A_R  = DSQRT(GM*P_R/R_R)
    
    
    
    
    
    
    
    
    
    !Part 13: 
    AC_BAR = 0.5*(A_L+A_R)
    
    !Part 14:
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
    !Part 15:
    M_L = UC_L/AC_BAR
    M_R = UC_R/AC_BAR
    M_BAR = 0.5*(M_L+M_R)
        
    !Part 16:
	IF(DABS(M_L)>1.0)THEN
	 M_PLUS = 0.5*(M_L+DABS(M_L))
	 
	ELSE
     M_PLUS = 0.25*(M_L+1.)*(M_L+1.)
	 
    END IF

    !Part 17:
	IF(DABS(M_R)>1.0)THEN
	 M_MINUS = 0.5*(M_R-DABS(M_R))
	 
	ELSE
     M_MINUS =-0.25*(M_R-1.)*(M_R-1.)
	 
    END IF

    !Part 18:
    MM = M_PLUS+M_MINUS
    
       
   !PART 19:
	IF(MM<=-1.0)THEN
        
        F1 = UC_R*R_R
        F2 = UC_R*R_R*U_R+P_R*NXX
        F3 = UC_R*R_R*V_R+P_R*NYY
        F4 = UC_R*(R_R*E_R+P_R)
        
    !Part 20:    
    ELSE IF(MM>=1.0)THEN
        
        F1 = UC_L*R_L
        F2 = UC_L*R_L*U_L+P_L*NXX
        F3 = UC_L*R_L*V_L+P_L*NYY
        F4 = UC_L*(R_L*E_L+P_L)
    
    !Part 21:
    ELSE
        
        !Part 22:
        !!ZHA-CUSP
        ALF_L = 2.0*((P_L/R_L)/((P_L/R_L)+(P_R/R_R)))
        ALF_R = 2.0*((P_R/R_R)/((P_L/R_L)+(P_R/R_R)))
        
        !Part 23:
        UP_L = AC_BAR*( 0.5*(M_L+DABS(M_L)) + ALF_L*(0.25*((M_L+1.0)**2)-0.5*(M_L+DABS(M_L))) )
        UM_R = AC_BAR*( 0.5*(M_R-DABS(M_R)) + ALF_R*(-0.25*((M_R-1.0)**2)-0.5*(M_R-DABS(M_R))) )
        
        !Part 24:
        RU_HALF = R_L*UP_L+R_R*UM_R
        
        !Part 25:
        F1_C = 0.5*(RU_HALF*(1.0+1.0)-DABS(RU_HALF)*(1.0-1.0))
        F2_C = 0.5*(RU_HALF*(U_L+U_R)-DABS(RU_HALF)*(U_R-U_L))
        F3_C = 0.5*(RU_HALF*(V_L+V_R)-DABS(RU_HALF)*(V_R-V_L))
        
        !Part 26:
        IZHA = 1
        
        !Part 27:
        IF (IZHA==1) THEN
            
            !!ZHA-CUSP2
            ALF_L = 2.0*((H_L/R_L)/((H_L/R_L)+(H_R/R_R)))
            ALF_R = 2.0*((H_R/R_R)/((H_L/R_L)+(H_R/R_R)))
            
            UP_L = AC_BAR*( 0.5*(M_L+DABS(M_L)) + ALF_L*(0.25*((M_L+1.0)**2)-0.5*(M_L+DABS(M_L))) )
            UM_R = AC_BAR*( 0.5*(M_R-DABS(M_R)) + ALF_R*(-0.25*((M_R-1.0)**2)-0.5*(M_R-DABS(M_R))) )
            
            RU_HALF = R_L*UP_L+R_R*UM_R
            
        END IF
        
        !Part 28:
        F4_C = 0.5*(RU_HALF*(E_L+E_R)-DABS(RU_HALF)*(E_R-E_L))
        
        !Part 29:
        ALF = 3.0/16.0
        P_PL = 0.25*((M_L+1.0)**2)*(2.0-M_L)+ALF*M_L*((M_L**2-1.0)**2)
        P_MR = 0.25*((M_R-1.0)**2)*(2.0+M_R)-ALF*M_R*((M_R**2-1.0)**2)
        
        !Part 30:
        F1_P = 0.0
        F2_P = P_PL*P_L*NXX+P_MR*P_R*NXX
        F3_P = P_PL*P_L*NYY+P_MR*P_R*NYY
        F4_P = 0.5*(P_L*(UC_L+AC_BAR)+P_R*(UC_R-AC_BAR))
        
        !Part 31:
        F1 = F1_C + F1_P 
        F2 = F2_C + F2_P 
        F3 = F3_C + F3_P 
        F4 = F4_C + F4_P 
                
	ENDIF

    
   !Part 32: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
                            
   !Part 33:                
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations Using ROE METHOD//!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: MARCH, 05, 2016                                                                //!
!// Developed by: S. KAVOUSI, Iran, FARS, MAMASANI                                       //!
!// Doc ID: MC2F002F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine HLL_METHOD_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L,P_L,E_L
 Real(8)::R_R,U_R,V_R,H_R,P_R,E_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,DEL_UC,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 
 
 Real(8):: R_BAR,A_BAR,P_BAR,P_PVRS,P_STAR,Q_L,Q_R,S_L,S_R,S_STAR,UL_COEF,UR_COEF
 Real(8):: U1_SL,U2_SL,U3_SL,U4_SL,U1_SR,U2_SR,U3_SR,U4_SR
 Real(8):: F1_SL,F2_SL,F3_SL,F4_SL,F1_SR,F2_SR,F3_SR,F4_SR
 Real(8):: F1_HLLC,F2_HLLC,F3_HLLC,F4_HLLC
 
 Real(8),DIMENSION(1:4)::DEL_W
 Real(8):: G1_HLLC,G2_HLLC,G3_HLLC,G4_HLLC
 Real(8)::G1_L,G2_L,G3_L,G4_L
 Real(8)::G1_R,G2_R,G3_R,G4_R
 
 Real(8)::P_SL,P_SR,ALP_L,BET_L,OMG_L,ALP_R,BET_R,OMG_R,Z_TR
 integer::ISTR
 
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,RHO_R
 Real(8)::PR,PHI,DX,DY
 
 
 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)
       R_L = Rho_L
       H_L = (R_L*E_L+P_L)/R_L
       A_L  = DSQRT(GM*P_L/R_L)

       
      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
       R_R = Rho_R
       H_R = (R_R*E_R+P_R)/R_R
       A_R  = DSQRT(GM*P_R/R_R)
    
    
   !Part 13: 
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
   !Part 14:
    ISTR = 5
   
   !Part 15:
    IF (ISTR==1) THEN
        
        S_L = UC_L-A_L
        S_R = UC_R+A_R
   !Part 16:     
    ELSE IF (ISTR==2) THEN
        
        S_L = DMIN1(UC_L-A_L,UC_R-A_R)
        S_R = DMAX1(UC_L+A_L,UC_R+A_R)
   !Part 17:     
    ELSE IF (ISTR==3) THEN
        
        R_RT = DSQRT(R_R/R_L)
        R_ROE = DSQRT(R_L*R_R)
        U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
        V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
        H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
        A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
        
        UC_ROE = U_ROE*NXX+V_ROE*NYY
        
        S_L = UC_ROE-A_ROE
        S_R = UC_ROE+A_ROE
   !Part 18:     
    ELSE IF (ISTR==4) THEN
        
        R_RT = DSQRT(R_R/R_L)
        R_ROE = DSQRT(R_L*R_R)
        U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
        V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
        H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
        A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
        UC_ROE = U_ROE*NXX+V_ROE*NYY
        
        S_L = DMIN1(UC_L-A_L,UC_ROE-A_ROE)
        S_R = DMAX1(UC_R+A_R,UC_ROE+A_ROE)
        
   !Part 19:     
    ELSE IF (ISTR==5) THEN
        
        R_BAR = 0.5*(R_L+R_R)                                                                        
        A_BAR = 0.5*(A_L+A_R)
        P_BAR = 0.5*(P_L+P_R)
        
        P_PVRS = P_BAR-0.5*(UC_R-UC_L)*R_BAR*A_BAR
        P_STAR = DMAX1(0.D0,P_PVRS)
        
        IF (P_STAR<=P_L) THEN
            
            Q_L = 1.0
            
        ELSE
            
            Q_L = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_L)-1.0))
            
        END IF
        
        IF (P_STAR<=P_R) THEN
            
            Q_R = 1.0
            
        ELSE
            
            Q_R = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_R)-1.0))
            
        END IF
        
        S_L = UC_L-A_L*Q_L
        S_R = UC_R+A_R*Q_R
        
        
   !Part 20:     
    ELSE IF (ISTR==6) THEN
        
        Z_TR = (GM-1.0)/(2.0*GM)
        
        P_STAR = ((A_L+A_R-0.5*(GM-1.0)*(U_R-U_L))/((A_L/(P_L**Z_TR))+(A_R/(P_R**Z_TR))))**(1.0/Z_TR)
        
    END IF
   
   !Part 21: 
    DEL_W(1) = WNP1(1,R)-WNP1(1,L)
    DEL_W(2) = WNP1(2,R)-WNP1(2,L)
    DEL_W(3) = WNP1(3,R)-WNP1(3,L)
    DEL_W(4) = WNP1(4,R)-WNP1(4,L)
    
   !Part 22: 
    IF (S_L>=0.0) THEN
        
        F1 = (UC_L*R_L               )
        F2 = (UC_L*R_L*U_L+P_L*NXX)
        F3 = (UC_L*R_L*V_L+P_L*NYY)
        F4 = (UC_L*(R_L*E_L+P_L)  )
   !Part 23:     
    ELSE IF (S_L<=0.0 .AND. S_R>=0.0) THEN
        
        F1_L = (UC_L*R_L               )
        F2_L = (UC_L*R_L*U_L+P_L*NXX)
        F3_L = (UC_L*R_L*V_L+P_L*NYY)
        F4_L = (UC_L*(R_L*E_L+P_L)  )
        
        F1_R = (UC_R*R_R               )
        F2_R = (UC_R*R_R*U_R+P_R*NXX)
        F3_R = (UC_R*R_R*V_R+P_R*NYY)
        F4_R = (UC_R*(R_R*E_R+P_R)  )
        
        F1 = ( (S_R*F1_L-S_L*F1_R+S_L*S_R*(DEL_W(1)))/(S_R-S_L) )
        F2 = ( (S_R*F2_L-S_L*F2_R+S_L*S_R*(DEL_W(2)))/(S_R-S_L) )
        F3 = ( (S_R*F3_L-S_L*F3_R+S_L*S_R*(DEL_W(3)))/(S_R-S_L) )
        F4 = ( (S_R*F4_L-S_L*F4_R+S_L*S_R*(DEL_W(4)))/(S_R-S_L) )
   !Part 24:     
    ELSE IF (S_R<=0) THEN
        
        F1 = (UC_R*R_R               )
        F2 = (UC_R*R_R*U_R+P_R*NXX)
        F3 = (UC_R*R_R*V_R+P_R*NYY)
        F4 = (UC_R*(R_R*E_R+P_R)  )
        
    END IF
    
     
   !Part 25: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
                            
   !Part 26:                
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations Using ROE METHOD//!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: MARCH, 05, 2016                                                                //!
!// Developed by: S. K                                                                   //!
!// Doc ID: MC2F002F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine HLLC_METHOD_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L,P_L,E_L
 Real(8)::R_R,U_R,V_R,H_R,P_R,E_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,DEL_UC,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 
 
 Real(8):: R_BAR,A_BAR,P_BAR,P_PVRS,P_STAR,Q_L,Q_R,S_L,S_R,S_STAR,UL_COEF,UR_COEF
 Real(8):: U1_SL,U2_SL,U3_SL,U4_SL,U1_SR,U2_SR,U3_SR,U4_SR
 Real(8):: F1_SL,F2_SL,F3_SL,F4_SL,F1_SR,F2_SR,F3_SR,F4_SR
 Real(8):: F1_HLLC,F2_HLLC,F3_HLLC,F4_HLLC
 
 Real(8),DIMENSION(1:4)::DEL_W
 Real(8):: G1_HLLC,G2_HLLC,G3_HLLC,G4_HLLC
 Real(8)::G1_L,G2_L,G3_L,G4_L
 Real(8)::G1_R,G2_R,G3_R,G4_R
 
 Real(8)::P_SL,P_SR,ALP_L,BET_L,OMG_L,ALP_R,BET_R,OMG_R,Z_TR,P_S
 integer::ISTR
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,RHO_R
 Real(8)::PR,PHI,DX,DY
 
 
 
 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)
       R_L = Rho_L
       H_L = (R_L*E_L+P_L)/R_L
       A_L  = DSQRT(GM*P_L/R_L)

       
      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
       R_R = Rho_R
       H_R = (R_R*E_R+P_R)/R_R
       A_R  = DSQRT(GM*P_R/R_R)
    
    
    
    
    
    
   !Part 13: 
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
   !Part 14:
    ISTR = 5
   
   !Part 15: 
    IF (ISTR==1) THEN
        
        S_L = UC_L-A_L
        S_R = UC_R+A_R
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
   !Part 16:     
    ELSE IF (ISTR==2) THEN
        
        S_L = DMIN1(UC_L-A_L,UC_R-A_R)
        S_R = DMAX1(UC_L+A_L,UC_R+A_R)
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
   
   !Part 17:     
    ELSE IF (ISTR==3) THEN
        
        R_RT = DSQRT(R_R/R_L)
        R_ROE = DSQRT(R_L*R_R)
        U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
        V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
        H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
        A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
        
        UC_ROE = U_ROE*NXX+V_ROE*NYY
        
        S_L = UC_ROE-A_ROE
        S_R = UC_ROE+A_ROE
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
   !Part 18:     
    ELSE IF (ISTR==4) THEN
        
        R_RT = DSQRT(R_R/R_L)
        R_ROE = DSQRT(R_L*R_R)
        U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
        V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
        H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
        A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
        
        UC_ROE = U_ROE*NXX+V_ROE*NYY
        
        S_L = DMIN1(UC_L-A_L,UC_ROE-A_ROE)
        S_R = DMAX1(UC_R+A_R,UC_ROE+A_ROE)
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
   !Part 19:     
    ELSE IF (ISTR==5) THEN
        
        R_BAR = 0.5*(R_L+R_R)                                                                        
        A_BAR = 0.5*(A_L+A_R)
        P_BAR = 0.5*(P_L+P_R)
        
        P_PVRS = P_BAR-0.5*(UC_R-UC_L)*R_BAR*A_BAR
        P_STAR = DMAX1(0.D0,P_PVRS)
        
        IF (P_STAR<=P_L) THEN
            
            Q_L = 1.0
            
        ELSE
            
            Q_L = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_L)-1.0))
            
        END IF
        
        IF (P_STAR<=P_R) THEN
            
            Q_R = 1.0
            
        ELSE
            
            Q_R = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_R)-1.0))
            
        END IF
        
        S_L = UC_L-A_L*Q_L
        S_R = UC_R+A_R*Q_R
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
    !Part 20:     
    ELSE IF (ISTR==6) THEN
        
        Z_TR = (GM-1.0)/(2.0*GM)
        
        P_STAR = ((A_L+A_R-0.5*(GM-1.0)*(U_R-U_L))/((A_L/(P_L**Z_TR))+(A_R/(P_R**Z_TR))))**(1.0/Z_TR)
       
        IF (P_STAR<=P_L) THEN
            
            Q_L = 1.0
            
        ELSE
            
            Q_L = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_L)-1.0))
            
        END IF
        
        IF (P_STAR<=P_R) THEN
            
            Q_R = 1.0
            
        ELSE
            
            Q_R = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_R)-1.0))
            
        END IF
        
        S_L = UC_L-A_L*Q_L
        S_R = UC_R+A_R*Q_R
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
    END IF
     
    !Part 21:     
     P_SL = P_L+R_L*(UC_L-S_L)*(UC_L-S_STAR)
     P_SR = P_R+R_R*(UC_R-S_R)*(UC_R-S_STAR)
     
     P_S = 0.5*(P_SL+P_SR)
     
    !Part 22:
     ALP_L = (S_L-UC_L)/(S_L-S_STAR)
     BET_L = ALP_L*(S_STAR-UC_L)
    
     ALP_R = (S_R-UC_R)/(S_R-S_STAR)
     BET_R = ALP_R*(S_STAR-UC_R)
     
    !Part 23: 
     U1_SL = ALP_L*(R_L) + (0.0)
     U2_SL = ALP_L*(R_L*U_L) + (R_L*BET_L*NXX)
     U3_SL = ALP_L*(R_L*V_L) + (R_L*BET_L*NYY)
     U4_SL = ALP_L*(R_L*E_L) + ((P_S*S_STAR-P_L*UC_L)/(S_L-S_STAR))
     
     U1_SR = ALP_R*(R_R) + (0.0)
     U2_SR = ALP_R*(R_R*U_R) + (R_R*BET_R*NXX)
     U3_SR = ALP_R*(R_R*V_R) + (R_R*BET_R*NYY)
     U4_SR = ALP_R*(R_R*E_R) + ((P_S*S_STAR-P_R*UC_R)/(S_R-S_STAR))
     
    !Part 24: 
     F1_L = (UC_L*R_L               )
     F2_L = (UC_L*R_L*U_L+P_L*NXX)
     F3_L = (UC_L*R_L*V_L+P_L*NYY)
     F4_L = (UC_L*(R_L*E_L+P_L)  )
    
     F1_R = (UC_R*R_R               )
     F2_R = (UC_R*R_R*U_R+P_R*NXX)
     F3_R = (UC_R*R_R*V_R+P_R*NYY)
     F4_R = (UC_R*(R_R*E_R+P_R)  )
     
    !Part 25: 
     IF (S_L>=0.0) THEN
         
         F1 = F1_L
         F2 = F2_L
         F3 = F3_L
         F4 = F4_L
         
    !Part 26:     
     ELSE IF (S_L<=0.0 .AND. S_STAR>=0.0) THEN
         
         F1 = F1_L+S_L*(U1_SL-R_L    )
         F2 = F2_L+S_L*(U2_SL-R_L*U_L)
         F3 = F3_L+S_L*(U3_SL-R_L*V_L)
         F4 = F4_L+S_L*(U4_SL-R_L*E_L)
         
    !Part 27:     
     ELSE IF (S_STAR<=0.0 .AND. S_R>=0.0) THEN
         
         F1 = F1_R+S_R*(U1_SR-R_R    )
         F2 = F2_R+S_R*(U2_SR-R_R*U_R)
         F3 = F3_R+S_R*(U3_SR-R_R*V_R)
         F4 = F4_R+S_R*(U4_SR-R_R*E_R)
         
    !Part 28:     
     ELSE IF (S_R<=0) THEN
         
         F1 = F1_R
         F2 = F2_R
         F3 = F3_R
         F4 = F4_R
         
     END IF
    
   !Part 29: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
                            
   !Part 30:                
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations Using ROE METHOD//!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: MARCH, 05, 2016                                                                //!
!// Developed by: N. msnkre, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F002F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ROE_EC_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*********************************************************************************************
  
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L,P_L
 Real(8)::R_R,U_R,V_R,H_R,P_R

 Real(8)::UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 
 Real(8)::Z1_BAR,Z2_BAR,Z3_BAR,Z4_BAR,Z1_LN,Z4_LN
 Real(8)::R_HAT,U_HAT,V_HAT,P1_HAT,P2_HAT,A_HAT,H_HAT,UC_HAT
 Real(8)::F1_CNV,F2_CNV,F3_CNV,F4_CNV
 Real(8),DIMENSION(1:4,1:4)::R_REC,RT_REC,JCB_REC
 Real(8),DIMENSION(1:4)::EV_REC,S_REC,DV_REC,DSP_REC
 INTEGER::IREC,JREC,KREC
 
 Real(8)::LOG_AVE_HO
 Real(8)::C1,C2
 
 INTEGER::IES
 Real(8)::DM_MAX,BET_EC2,ALF_MIN,ALF_MAX,ALF_EC2,ALF_EC1
 
 
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,E_L,RHO_R,E_R
 Real(8)::PR,PHI,DX,DY

 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do
 
 

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)
       R_L = Rho_L
       H_L = (R_L*E_L+P_L)/R_L
       A_L  = DSQRT(GM*P_L/R_L)

       
      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
       R_R = Rho_R
       H_R = (R_R*E_R+P_R)/R_R
       A_R  = DSQRT(GM*P_R/R_R)
    
    
    
   !Part 13:
    Z1_BAR = 0.5*(DSQRT(R_L/P_L)+DSQRT(R_R/P_R))
    Z2_BAR = 0.5*(DSQRT(R_L/P_L)*U_L+DSQRT(R_R/P_R)*U_R)
    Z3_BAR = 0.5*(DSQRT(R_L/P_L)*V_L+DSQRT(R_R/P_R)*V_R)              
    Z4_BAR = 0.5*(DSQRT(R_L*P_L)+DSQRT(R_R*P_R))                  
   
   !Part 14: 
    C1 = DSQRT(R_L/P_L)
    C2 = DSQRT(R_R/P_R)
    Z1_LN = LOG_AVE_HO(C1,C2)
    
    C1 = DSQRT(R_L*P_L)
    C2 = DSQRT(R_R*P_R)
    Z4_LN = LOG_AVE_HO(C1,C2)
    
   !Part 15: 
    R_HAT = Z1_BAR*Z4_LN
    U_HAT = Z2_BAR/Z1_BAR
    V_HAT = Z3_BAR/Z1_BAR
    P1_HAT = Z4_BAR/Z1_BAR
    P2_HAT = ((GM+1)*(Z4_LN))/(2.0*GM*Z1_LN)+((GM-1)*(Z4_BAR))/(2.0*GM*Z1_BAR)
    
    A_HAT = DSQRT((GM*P2_HAT)/R_HAT)
    H_HAT = (A_HAT**2)/(GM-1.0)+0.5*(U_HAT**2+V_HAT**2)
    
    UC_HAT = U_HAT*NXX+V_HAT*NYY

   !Part 16:  
    F1_CNV = UC_HAT*R_HAT
    F2_CNV = UC_HAT*R_HAT*U_HAT+P1_HAT*NXX
    F3_CNV = UC_HAT*R_HAT*V_HAT+P1_HAT*NYY
    F4_CNV = UC_HAT*R_HAT*H_HAT
    
   !Part 17:  
    R_REC(1,1) = 1.0
    R_REC(1,2) = 1.0
    R_REC(1,3) = 0.0
    R_REC(1,4) = 1.0
    
    R_REC(2,1) = U_HAT-A_HAT*NXX
    R_REC(2,2) = U_HAT
    R_REC(2,3) = NYY
    R_REC(2,4) = U_HAT+A_HAT*NXX
    
    R_REC(3,1) = V_HAT-A_HAT*NYY
    R_REC(3,2) = V_HAT
    R_REC(3,3) =-NXX
    R_REC(3,4) = V_HAT+A_HAT*NYY
    
    R_REC(4,1) = H_HAT-A_HAT*UC_HAT
    R_REC(4,2) = 0.5*(U_HAT**2+V_HAT**2)
    R_REC(4,3) = U_HAT*NYY-V_HAT*NXX
    R_REC(4,4) = H_HAT+A_HAT*UC_HAT
    
   !Part 18:  
    IES=2
    IF (IES==0) THEN
        
        !Part 19: 
        !!ENTROPY STABLE SCHEME, ROE-ES
        EV_REC(1) = DABS(UC_HAT-A_HAT)
        EV_REC(2) = DABS(UC_HAT)
        EV_REC(3) = DABS(UC_HAT)  
        EV_REC(4) = DABS(UC_HAT+A_HAT)
        
    ELSE IF (IES==1) THEN
        
        !Part 20: 
        !!ENTROPY CONSISTENT#1 SCHEME, ROE-EC1
        UC_L = U_L*NXX+V_L*NYY
        UC_R = U_R*NXX+V_R*NYY
        
        !Part 21:
        ALF_EC1 = 1.0/6.0       
        
        EV_REC(1) = DABS(UC_HAT-A_HAT)+(ALF_EC1)*DABS( (UC_R-A_R)-(UC_L-A_L) )
        EV_REC(2) = DABS(UC_HAT)
        EV_REC(3) = DABS(UC_HAT)  
        EV_REC(4) = DABS(UC_HAT+A_HAT)+(ALF_EC1)*DABS( (UC_R+A_R)-(UC_L+A_L) )
        
    ELSE IF (IES==2) THEN
        
        !Part 22: 
        !!ENTROPY CONSISTENT#2 SCHEME, ROE-EC2
        ALF_MIN = 1.0/6.0
        ALF_MAX = 2.0
        DM_MAX = 0.5
        BET_EC2 = 1.0/6.0
        
        !Part 23:
        UC_L = U_L*NXX+V_L*NYY
        UC_R = U_R*NXX+V_R*NYY
        
        !Part 24:
        M_L = DSQRT( (U_L**2+V_L**2)/(GM*P_L/R_L) )
        M_R = DSQRT( (U_R**2+V_R**2)/(GM*P_R/R_R) ) 
        
        !Part 25: 
        ALF_EC2 = (ALF_MAX-ALF_MIN)*( DMAX1(0.0,DSIGN(DM_MAX,(M_R-M_L))) ) + ALF_MIN               
        
        EV_REC(1) = (1.0+BET_EC2)*DABS(UC_HAT-A_HAT)+(ALF_EC2)*DABS( (UC_R-A_R)-(UC_L-A_L) )
        EV_REC(2) = DABS(UC_HAT)
        EV_REC(3) = DABS(UC_HAT)  
        EV_REC(4) = (1.0+BET_EC2)*DABS(UC_HAT+A_HAT)+(ALF_EC2)*DABS( (UC_R+A_R)-(UC_L+A_L) )
        
        
    END IF
    
    !Part 26:
    S_REC(1) = (R_HAT)/(2.0*GM)
    S_REC(2) = ((GM-1)*R_HAT)/GM
    S_REC(3) = P1_HAT                              
    S_REC(4) = (R_HAT)/(2.0*GM)
    
    !Part 27:
    DV_REC(1) = ( (GM-LOG(P_R/(R_R**GM)))/(GM-1.0)-0.5*(R_R/P_R)*(U_R**2+V_R**2) ) - ( (GM-LOG(P_L/(R_L**GM)))/(GM-1.0)-0.5*(R_L/P_L)*(U_L**2+V_L**2) )
    DV_REC(2) = ( (R_R/P_R)*U_R ) - ( (R_L/P_L)*U_L )
    DV_REC(3) = ( (R_R/P_R)*V_R ) - ( (R_L/P_L)*V_L )
    DV_REC(4) = -(R_R/P_R) + (R_L/P_L)
    
    !Part 28:
    RT_REC(1,1) = 1.0                                                 
    RT_REC(1,2) = U_HAT-A_HAT*NXX
    RT_REC(1,3) = V_HAT-A_HAT*NYY
    RT_REC(1,4) = H_HAT-A_HAT*UC_HAT
    
    RT_REC(2,1) = 1.0
    RT_REC(2,2) = U_HAT
    RT_REC(2,3) = V_HAT
    RT_REC(2,4) = 0.5*(U_HAT**2+V_HAT**2)
    
    RT_REC(3,1) = 0.0
    RT_REC(3,2) = NYY
    RT_REC(3,3) =-NXX
    RT_REC(3,4) = U_HAT*NYY-V_HAT*NXX
    
    RT_REC(4,1) = 1.0
    RT_REC(4,2) = U_HAT+A_HAT*NXX
    RT_REC(4,3) = V_HAT+A_HAT*NYY
    RT_REC(4,4) = H_HAT+A_HAT*UC_HAT
    
    !Part 29:
    JCB_REC(:,:) = 0
    
    !Part 30:
    DO IREC=1,4
        DO JREC=1,4
            DO KREC=1,4
                
                JCB_REC(IREC,JREC) = JCB_REC(IREC,JREC)+R_REC(IREC,KREC)*EV_REC(KREC)*S_REC(KREC)*RT_REC(KREC,JREC)         !!!!!!
                
            END DO
        END DO
    END DO
            
   !Part 31: 
    DSP_REC(:) = 0
    
    !Part 32:
    DO IREC=1,4
        DO JREC=1,4
            
            DSP_REC(IREC) = DSP_REC(IREC)+JCB_REC(IREC,JREC)*DV_REC(JREC)
            
        END DO
    END DO
   
   !Part 33: 
    F1 = F1_CNV-0.5*DSP_REC(1)
    F2 = F2_CNV-0.5*DSP_REC(2)
    F3 = F3_CNV-0.5*DSP_REC(3)
    F4 = F4_CNV-0.5*DSP_REC(4)
    
   !Part 34: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
    
   !Part 35: 
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
End SUBROUTINE
    
    
!Part 36:    
Real(8) FUNCTION LOG_AVE_HO(PRM_L,PRM_R)
   
   IMPLICIT NONE
   INTENT(IN   )::PRM_L,PRM_R
   Real(8)::PRM_L,PRM_R
   Real(8)::EPS,ZETA,F0,U0,F1
   
   EPS = 0.01
   
   ZETA = PRM_L/PRM_R
   
   F0 = (ZETA-1.0)/(ZETA+1.0)
   
   U0 = F0**2
   
   IF (U0<EPS) THEN
       
       F1 = 1+(U0/3)+((U0**2)/5)+((U0**3)/7)
       
   ELSE
       
       F1 = LOG(ZETA)/(2.0*F0)
       
   END IF
   
   LOG_AVE_HO = (PRM_L+PRM_R)/(2.0*F1)
   
    
   RETURN
   
END FUNCTION LOG_AVE_HO
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations Using ROE METHOD//!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: MARCH, 05, 2016                                                                //!
!// Developed by: S. KAVOUSI, Iran, FARS, MAMASANI                                       //!
!// Doc ID: MC2F002F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ROE_METHOD_ENT_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L
 Real(8)::R_R,U_R,V_R,H_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,DEL_UC,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 Real(8)::DEL_R,DEL_U,DEL_V,DEL_P
 Real(8)::ROE_COEF1,ROE_COEF2,ROE_COEF3,ROE_COEF4,ROE_COEF5
 Real(8)::F1_ROE1,F2_ROE1,F3_ROE1,F4_ROE1
 Real(8)::F1_ROE4,F2_ROE4,F3_ROE4,F4_ROE4
 Real(8)::F1_ROE5,F2_ROE5,F3_ROE5,F4_ROE5
 
 INTEGER::IENT
 Real(8)::EV1_R,EV2_R,EV3_R,EV1_L,EV2_L,EV3_L,EV1_ROE,EV2_ROE,EV3_ROE
 Real(8)::EPS1,EPS2,EPS3,EV1_NEW,EV2_NEW,EV3_NEW
 
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,P_L,E_L,RHO_R,P_R,E_R
 Real(8)::PR,PHI,DX,DY



 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2
     
     
     
   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA
    
    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)
       R_L = Rho_L
       H_L = (R_L*E_L+P_L)/R_L
       A_L  = DSQRT(GM*P_L/R_L)

       
      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
       R_R = Rho_R
       H_R = (R_R*E_R+P_R)/R_R
       A_R  = DSQRT(GM*P_R/R_R)
    
    
    
    
    
   !Part 13: 
    R_RT = DSQRT(R_R/R_L)

    R_ROE = DSQRT(R_L*R_R)
    U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
    V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
    H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
    
    A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
    
   !Part 14: 
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
    UC_ROE = U_ROE*NXX+V_ROE*NYY
    
    DEL_UC = (U_R-U_L)*NXX+(V_R-V_L)*NYY
    
   !Part 15: 
    F1_L = UC_L*R_L
    F2_L = UC_L*R_L*U_L+P_L*NXX
    F3_L = UC_L*R_L*V_L+P_L*NYY
    F4_L = UC_L*(R_L*E_L+P_L)
    
    F1_R = UC_R*R_R
    F2_R = UC_R*R_R*U_R+P_R*NXX
    F3_R = UC_R*R_R*V_R+P_R*NYY
    F4_R = UC_R*(R_R*E_R+P_R)
    
   !Part 16: 
    DEL_R = R_R-R_L
    DEL_U = U_R-U_L
    DEL_V = V_R-V_L
    DEL_P = P_R-P_L
    
    
   !Part 17: 
    IENT = 2
    !!ENTROPY FIX CORRECTIONS!
    IF (IENT /= 0) THEN
        
        !Part 18:
        EV1_R = UC_R
        EV2_R = UC_R+A_R*DSQRT(NXX*NXX+NYY*NYY)
        EV3_R = UC_R-A_R*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_L = UC_L
        EV2_L = UC_L+A_L*DSQRT(NXX*NXX+NYY*NYY)
        EV3_L = UC_L-A_L*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_ROE = UC_ROE
        EV2_ROE = UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        EV3_ROE = UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        
        !Part 19:
        IF (IENT == 1) THEN
            
            !!ENTROPY FIX KERMANI#1!
            EPS1 = 2.D0*DMAX1(0.D0,(EV1_R-EV1_L))
            EPS2 = 2.D0*DMAX1(0.D0,(EV2_R-EV2_L))
            EPS3 = 2.D0*DMAX1(0.D0,(EV3_R-EV3_L))
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 20:    
        ELSE IF (IENT == 2) THEN
            
            !!ENTROPY FIX KERMANI#2!
            EPS1 = 4.D0*DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = 4.D0*DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = 4.D0*DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 21:    
        ELSE IF (IENT == 3) THEN
            
            !!ENTROPY FIX HARTEN-HYMAN#1!
            EPS1 = DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 22:    
        ELSE IF (IENT == 4) THEN
            
            !!ENTROPY FIX HARTEN-HYMAN#2!
            EPS1 = DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = EPS1
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = EPS2
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = EPS3
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 23:
        ELSE IF (IENT == 5) THEN
            
            !!ENTROPY FIX HOFFMAN-CHIANG!
            EPS1 = 0.05          !EPS SHOULD BE: 0<EPS<0.125
            EPS2 = 0.05           !EPS SHOULD BE: 0<EPS<0.125
            EPS3 = 0.05           !EPS SHOULD BE: 0<EPS<0.125
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = EPS1
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = EPS2
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = EPS3
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
        END IF
        
        
    ELSE
        !Part 24:
        EV1_NEW = UC_ROE
        EV2_NEW = UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY)  
        EV3_NEW = UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY)  
    
    END IF
    
    
   !Part 25: 
    ROE_COEF1 = DABS(EV1_NEW)*(DEL_R-DEL_P/(A_ROE*A_ROE))
    ROE_COEF2 = DABS(EV1_NEW)*R_ROE
    ROE_COEF3 = 0.5*(U_ROE*U_ROE+V_ROE*V_ROE)
    ROE_COEF4 = DABS(EV2_NEW)*((DEL_P+(R_ROE*A_ROE*DEL_UC))/(2.0*A_ROE*A_ROE))
    ROE_COEF5 = DABS(EV3_NEW)*((DEL_P-(R_ROE*A_ROE*DEL_UC))/(2.0*A_ROE*A_ROE))
    
   !Part 26: 
    F1_ROE1 = ROE_COEF1*1.0+ROE_COEF2*0.0
    F2_ROE1 = ROE_COEF1*U_ROE+ROE_COEF2*(DEL_U-NXX*DEL_UC)
    F3_ROE1 = ROE_COEF1*V_ROE+ROE_COEF2*(DEL_V-NYY*DEL_UC)
    F4_ROE1 = ROE_COEF1*ROE_COEF3+ROE_COEF2*(U_ROE*DEL_U+V_ROE*DEL_V-UC_ROE*DEL_UC)
    
    F1_ROE4 = ROE_COEF4*(1.0)
    F2_ROE4 = ROE_COEF4*(U_ROE+A_ROE*NXX)
    F3_ROE4 = ROE_COEF4*(V_ROE+A_ROE*NYY)
    F4_ROE4 = ROE_COEF4*(H_ROE+A_ROE*UC_ROE)
    
    F1_ROE5 = ROE_COEF5*(1.0)
    F2_ROE5 = ROE_COEF5*(U_ROE-A_ROE*NXX)
    F3_ROE5 = ROE_COEF5*(V_ROE-A_ROE*NYY)
    F4_ROE5 = ROE_COEF5*(H_ROE-A_ROE*UC_ROE)
    
    
   !Part 27: 
    F1 = 0.5*(F1_L+F1_R-(F1_ROE1+F1_ROE4+F1_ROE5))
    F2 = 0.5*(F2_L+F2_R-(F2_ROE1+F2_ROE4+F2_ROE5))
    F3 = 0.5*(F3_L+F3_R-(F3_ROE1+F3_ROE4+F3_ROE5))
    F4 = 0.5*(F4_L+F4_R-(F4_ROE1+F4_ROE4+F4_ROE5))
    
   !Part 28: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
    
   !Part 29: 
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
    
    
    
    
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!////////////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Convection Terms of 2D Mean Flow Equations Using ROE_METHOD     //!
!//                                                                                            //!
!// Version: V1                                                                                //!
!// Date: October, 12, 2014                                                                    //!
!// Developed by: N. msnkre, Iran, Tehran, OpenFlows@chmail.ir                                 //!
!// Doc ID: MC2F008F1                                                                          //!
!//                                                                                            //!
!// The Program Is Available Through The Website: www.MarketCode.ir                            //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                       //!
!////////////////////////////////////////////////////////////////////////////////////////////////!
!*************************************************************************************************
 Subroutine ROE_METHOD_JCB_ENT_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con,GWNP1,Xc,Yc,Limit,X,Y)
 Implicit None
!*************************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Limit,X,Y
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 
 Real(8)::R_L,U_L,V_L,H_L
 Real(8)::R_R,U_R,V_R,H_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 Real(8)::ROE_COEF1,ROE_COEF2,ROE_COEF3
 
 INTEGER::IROE,JROE,KROE
 Real(8),DIMENSION(1:4)::DEL_W,ADW_ROE,EV_ROE
 Real(8),DIMENSION(1:4,1:4)::REV_ROE,LEV_ROE,JCB_ROE
 
 INTEGER::IENT
 Real(8)::EV1_R,EV2_R,EV3_R,EV1_L,EV2_L,EV3_L,EV1_ROE,EV2_ROE,EV3_ROE
 Real(8)::EPS1,EPS2,EPS3,EV1_NEW,EV2_NEW,EV3_NEW
 
 
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 
 Real(8)::RHO_L,P_L,E_L,RHO_R,P_R,E_R
 Real(8)::PR,PHI,DX,DY
 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA
    
    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    

      !Part 11:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro    
       Pr = P(L) 
    
       DX = XM_EDG - XC(L)
       DY = YM_EDG - YC(L)     

       Phi=Limit(1,L)
	   Call Recons2Ord(Dim,1,L,Ro,GWNP1,Phi,DX,DY,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord(Dim,2,L,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord(Dim,3,L,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,L)         
	   Call Recons2Ord(Dim,4,L,Pr,GWNP1,Phi,DX,DY,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       P_L   = Pr
       
       E_L = P_L/(Rho_L*(GM-1.0)) + 0.5*(U_L*U_L + V_L*V_L)
       R_L = Rho_L
       H_L = (R_L*E_L+P_L)/R_L
       A_L  = DSQRT(GM*P_L/R_L)

       
      !Part 12:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro    
       Pr = P(R) 
    
       DX = XM_EDG - XC(R)
       DY = YM_EDG - YC(R)     

       Phi=Limit(1,R)
	   Call Recons2Ord(Dim,1,R,Ro,GWNP1,Phi,DX,DY,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord(Dim,2,R,U,GWNP1,Phi,DX,DY,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord(Dim,3,R,V,GWNP1,Phi,DX,DY,V)
       Phi=Limit(4,R)         
	   Call Recons2Ord(Dim,4,R,Pr,GWNP1,Phi,DX,DY,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V  
       P_R   = Pr
       
       E_R = P_R/(Rho_R*(GM-1.0)) + 0.5*(U_R*U_R + V_R*V_R) 
       R_R = Rho_R
       H_R = (R_R*E_R+P_R)/R_R
       A_R  = DSQRT(GM*P_R/R_R)
    
    
    
    
    
    
    
    
    
    
    
    
   !Part 13:
    R_RT = DSQRT(R_R/R_L)

    R_ROE = DSQRT(R_L*R_R)
    U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
    V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
    H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
    
    A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
    
   !Part 14:
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY    
    
    UC_ROE = U_ROE*NXX+V_ROE*NYY
    
   !Part 15:    
    F1_L = UC_L*R_L
    F2_L = UC_L*R_L*U_L+P_L*NXX
    F3_L = UC_L*R_L*V_L+P_L*NYY
    F4_L = UC_L*(R_L*E_L+P_L)
    
    F1_R = UC_R*R_R
    F2_R = UC_R*R_R*U_R+P_R*NXX
    F3_R = UC_R*R_R*V_R+P_R*NYY
    F4_R = UC_R*(R_R*E_R+P_R)
   
    
   !Part 16:    
    ROE_COEF1 = R_ROE/(DSQRT(2.D0)*A_ROE)
    ROE_COEF2 = 1/(DSQRT(2.D0)*R_ROE*A_ROE)
    ROE_COEF3 = 0.5*(GM-1)*(U_ROE*U_ROE+V_ROE*V_ROE)
    
   !Part 17:    
    REV_ROE(1,1) = 1
    REV_ROE(1,2) = 0
    REV_ROE(1,3) = ROE_COEF1
    REV_ROE(1,4) = ROE_COEF1
    
    REV_ROE(2,1) = U_ROE
    REV_ROE(2,2) = R_ROE*NYY
    REV_ROE(2,3) = ROE_COEF1*(U_ROE+A_ROE*NXX)
    REV_ROE(2,4) = ROE_COEF1*(U_ROE-A_ROE*NXX)
    
    REV_ROE(3,1) = V_ROE
    REV_ROE(3,2) = -R_ROE*NXX
    REV_ROE(3,3) = ROE_COEF1*(V_ROE+A_ROE*NYY)
    REV_ROE(3,4) = ROE_COEF1*(V_ROE-A_ROE*NYY)
    
    REV_ROE(4,1) = (ROE_COEF3)/(GM-1)
    REV_ROE(4,2) = R_ROE*(U_ROE*NYY-V_ROE*NXX)
    REV_ROE(4,3) = ROE_COEF1*((ROE_COEF3+A_ROE*A_ROE)/(GM-1)+A_ROE*UC_ROE)
    REV_ROE(4,4) = ROE_COEF1*((ROE_COEF3+A_ROE*A_ROE)/(GM-1)-A_ROE*UC_ROE)
    
    
   !Part 18:    
    LEV_ROE(1,1) = 1-(ROE_COEF3/(A_ROE*A_ROE))
    LEV_ROE(1,2) = ((GM-1)*U_ROE)/(A_ROE*A_ROE)
    LEV_ROE(1,3) = ((GM-1)*V_ROE)/(A_ROE*A_ROE)
    LEV_ROE(1,4) = -(GM-1)/(A_ROE*A_ROE)
    
    LEV_ROE(2,1) = -(U_ROE*NYY-V_ROE*NXX)/R_ROE
    LEV_ROE(2,2) = NYY/R_ROE
    LEV_ROE(2,3) = -NXX/R_ROE
    LEV_ROE(2,4) = 0
    
    LEV_ROE(3,1) = ROE_COEF2*(ROE_COEF3-A_ROE*UC_ROE)
    LEV_ROE(3,2) = ROE_COEF2*(A_ROE*NXX-(GM-1)*U_ROE)
    LEV_ROE(3,3) = ROE_COEF2*(A_ROE*NYY-(GM-1)*V_ROE)
    LEV_ROE(3,4) = ROE_COEF2*(GM-1)
    
    LEV_ROE(4,1) = ROE_COEF2*(ROE_COEF3+A_ROE*UC_ROE)
    LEV_ROE(4,2) = -ROE_COEF2*(A_ROE*NXX+(GM-1)*U_ROE)
    LEV_ROE(4,3) = -ROE_COEF2*(A_ROE*NYY+(GM-1)*V_ROE)
    LEV_ROE(4,4) = ROE_COEF2*(GM-1)
    
    
    
   !Part 19:
    IENT = 2
    !!ENTROPY FIX CORRECTIONS!
    IF (IENT /= 0) THEN
        
        !Part 20:
        EV1_R = UC_R
        EV2_R = UC_R+A_R*DSQRT(NXX*NXX+NYY*NYY)
        EV3_R = UC_R-A_R*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_L = UC_L
        EV2_L = UC_L+A_L*DSQRT(NXX*NXX+NYY*NYY)
        EV3_L = UC_L-A_L*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_ROE = UC_ROE
        EV2_ROE = UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        EV3_ROE = UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        
        !Part 21:
        IF (IENT == 1) THEN
            
            !!ENTROPY FIX KERMANI#1!
            EPS1 = 2.D0*DMAX1(0.D0,(EV1_R-EV1_L))
            EPS2 = 2.D0*DMAX1(0.D0,(EV2_R-EV2_L))
            EPS3 = 2.D0*DMAX1(0.D0,(EV3_R-EV3_L))
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 22:
        ELSE IF (IENT == 2) THEN
            
            !!ENTROPY FIX KERMANI#2!
            EPS1 = 4.D0*DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = 4.D0*DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = 4.D0*DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 23:    
        ELSE IF (IENT == 3) THEN
            
            !!ENTROPY FIX HARTEN-HYMAN#1!
            EPS1 = DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 24:
        ELSE IF (IENT == 4) THEN
            
            !!ENTROPY FIX HARTEN-HYMAN#2!
            EPS1 = DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = EPS1
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = EPS2
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = EPS3
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 25:
        ELSE IF (IENT == 5) THEN
            
            !!ENTROPY FIX HOFFMAN-CHIANG!
            EPS1 = 0.1           !EPS SHOULD BE: 0<EPS<0.125
            EPS2 = 0.1           !EPS SHOULD BE: 0<EPS<0.125
            EPS3 = 0.1           !EPS SHOULD BE: 0<EPS<0.125
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = EPS1
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = EPS2
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = EPS3
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        
        END IF
        
        
    ELSE
        !Part 26:
        EV_ROE(1) = DABS(UC_ROE)
        EV_ROE(2) = DABS(UC_ROE)
        EV_ROE(3) = DABS(UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY))  
        EV_ROE(4) = DABS(UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY))  
    
    END IF
    
   
        
    
   !Part 27:    
    JCB_ROE(:,:) = 0
    
    !Part 28:
    DO IROE=1,4
        DO JROE=1,4
            DO KROE=1,4
                
                JCB_ROE(IROE,JROE) = JCB_ROE(IROE,JROE)+REV_ROE(IROE,KROE)*EV_ROE(KROE)*LEV_ROE(KROE,JROE)         !!!!!!
                
            END DO
        END DO
    END DO
            
    
   !Part 29:
    DEL_W(1) = R_R-R_L                   
    DEL_W(2) = R_R*U_R-R_L*U_L 
    DEL_W(3) = R_R*V_R-R_L*V_L 
    DEL_W(4) = R_R*E_R-R_L*E_L 
   
    
   !Part 30:    
    ADW_ROE(:) = 0
    
    !Part 31:
    DO IROE=1,4
        DO JROE=1,4
            
            ADW_ROE(IROE) = ADW_ROE(IROE)+JCB_ROE(IROE,JROE)*DEL_W(JROE)
            
        END DO
    END DO
    
   !Part 32:
    F1 = 0.5*(F1_L+F1_R-ADW_ROE(1))*DAA
    F2 = 0.5*(F2_L+F2_R-ADW_ROE(2))*DAA
    F3 = 0.5*(F3_L+F3_R-ADW_ROE(3))*DAA
    F4 = 0.5*(F4_L+F4_R-ADW_ROE(4))*DAA
    
   !Part 33:    
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4
    
   !Part 34:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
