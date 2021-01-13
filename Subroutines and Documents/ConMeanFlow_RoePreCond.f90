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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_RoePreCond(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Minf,Con)
 Implicit None
!*************************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Minf
 Intent(Out  )::Con

 Integer::Dim,J,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,V_tot,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mach,M,Mm,Minf,M2,Pm,Nxx,Nyy,DAA,kk,a
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 
 Real(8)::R_L,U_L,V_L,H_L
 Real(8)::R_R,U_R,V_R,H_R
 Real(8)::R_RT,d1,d2,d3,d4,d5,d6,d7,d8,d9
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,M_ROE,H_ROE
 Real(8)::UC_ROE,UC_L,UC_R,VC_ROE,VC_L,VC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 Real(8)::Alpha1,COE1,COE2
 Real(8)::T,X,phi
 
 Real(8),Dimension(1:4)::ADW_ROE,D
 Real(8),Dimension(1:4,1:4)::Tur
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
	
   !Part 11:
    R_L = WNP1(1,L)
    U_L = WNP1(2,L)/WNP1(1,L)
    V_L = WNP1(3,L)/WNP1(1,L)
    
    R_R = WNP1(1,R)
    U_R = WNP1(2,R)/WNP1(1,R)
    V_R = WNP1(3,R)/WNP1(1,R)
   
   !Part 12:
    H_L = (WNP1(4,L)+P(L))/WNP1(1,L)
    H_R = (WNP1(4,R)+P(R))/WNP1(1,R)
    
    A_L  = DSQRT(GM*P(L)/WNP1(1,L))
	A_R  = DSQRT(GM*P(R)/WNP1(1,R))
    
   !Part 13:
    R_RT = DSQRT(R_R/R_L)

    R_ROE = DSQRT(R_L*R_R)
    U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
    V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
    H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
    
    A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
	M_ROE=(U_ROE*U_ROE+V_ROE*V_ROE)/(A_ROE*A_ROE)

   !Part 14:
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY 
	
	VC_L = V_L*NXX-U_L*NYY
   	VC_R = V_R*NXX-U_R*NYY
    
    UC_ROE = U_ROE*NXX+V_ROE*NYY
	VC_ROE = V_ROE*NXX-U_ROE*NYY

   !Part 15:    
    F1_L = UC_L*R_L
    F2_L = UC_L*WNP1(2,L)+P(L)*NXX
    F3_L = UC_L*WNP1(3,L)+P(L)*NYY
    F4_L = UC_L*(WNP1(4,L)+P(L))
    
    F1_R = UC_R*R_R
    F2_R = UC_R*WNP1(2,R)+P(R)*NXX
    F3_R = UC_R*WNP1(3,R)+P(R)*NYY
    F4_R = UC_R*(WNP1(4,R)+P(R))

   !Part 16:
	M2=max(Minf*Minf,M_ROE*M_ROE)
    T=min(M2,1.0)

   !Part 17:
	X=DSQRT(((1-T)*UC_ROE)**2+4*T*A_ROE*A_ROE)
	COE1=0.5*((1+T)*UC_ROE+X)-UC_ROE*T
	COE2=0.5*((1+T)*UC_ROE-X)-UC_ROE*T
	Alpha1=WNP1(1,R)-WNP1(1,L)-(P(R)-P(L))/(A_ROE*A_ROE)
		

   !Part 18:
	d1=UC_ROE*(1+T)/X
	IF (UC_ROE/=0) then		
    	d2=UC_ROE*UC_ROE*(T-1+2*A_ROE*A_ROE/(UC_ROE*UC_ROE))/X
        d9=UC_ROE*UC_ROE*(-T+1+2*A_ROE*A_ROE*T/(UC_ROE*UC_ROE))/X
	else
	    d2=2*A_ROE*A_ROE/X
		d9=2*A_ROE*A_ROE*T/X
	end IF
	d3=U_ROE*d1+NXX*d9
    d4=U_ROE*d2-COE1*COE2*NXX*d1/T
	d5=V_ROE*d1+NYY*d9
	d6=V_ROE*d2-COE1*COE2*NYY*d1/T
	d7=H_ROE*d1+UC_ROE*d9
	d8=H_ROE*d2-COE1*COE2*UC_ROE*d1/T
		

   !Part 19:
	ADW_ROE(1) = ABS(UC_ROE)*Alpha1+d2*(P(R)-P(L))/(A_ROE*A_ROE)+d1*R_ROE*(UC_R-UC_L)
	ADW_ROE(2) = ABS(UC_ROE)*Alpha1*U_ROE-R_ROE*ABS(UC_ROE)*NYY*(VC_R-VC_L)+d4*(P(R)-P(L))/(A_ROE*A_ROE)+d3*R_ROE*(UC_R-UC_L)
	ADW_ROE(3) = ABS(UC_ROE)*Alpha1*V_ROE+R_ROE*ABS(UC_ROE)*NXX*(VC_R-VC_L)+d6*(P(R)-P(L))/(A_ROE*A_ROE)+d5*R_ROE*(UC_R-UC_L)
    ADW_ROE(4) = ABS(UC_ROE)*Alpha1*(V_ROE*V_ROE+U_ROE*U_ROE)/2+R_ROE*ABS(UC_ROE)*VC_ROE*(VC_R-VC_L)+d8*(P(R)-P(L))/(A_ROE*A_ROE)+d7*R_ROE*(UC_R-UC_L)
  	    
   
   !Part 20:
    F1 = 0.5*(F1_L+F1_R-ADW_ROE(1))*DAA
    F2 = 0.5*(F2_L+F2_R-ADW_ROE(2))*DAA
    F3 = 0.5*(F3_L+F3_R-ADW_ROE(3))*DAA
    F4 = 0.5*(F4_L+F4_R-ADW_ROE(4))*DAA
    
   !Part 21:    
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4
    
   !Part 22:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    
 End Do

!Part 23:
 Do J=1,NC
  
 !Part 24:
  U=WNP1(2,J)/WNP1(1,J)
  V=WNP1(3,J)/WNP1(1,J)
  V_tot=U*U+V*V
  Mach = DSqrt( (U*U+V*V) / (GM*P(J)/WNP1(1,J)) )
  a=DSQRT(GM*P(J)/WNP1(1,J))

 !Part 25:
  M=Mach*Mach
  m2=max(Minf*Minf,M)
  T=min(m2,1.0)
  phi=(GM-1)*(T-1)
  
 !Part 26:
  Tur(1,1)=phi*M/2+1
  Tur(1,2)=-phi*U/(a*a)
  Tur(1,3)=-phi*V/(a*a)
  Tur(1,4)=phi/(a*a)

  Tur(2,1)=phi*U*M/2
  Tur(2,2)=-phi*U*U/(a*a)+1
  Tur(2,3)=-phi*U*V/(a*a)
  Tur(2,4)=phi*U/(a*a)

  Tur(3,1)=phi*V*M/2
  Tur(3,2)=-phi*U*V/(a*a)
  Tur(3,3)=-phi*V*V/(a*a)+1
  Tur(3,4)=phi*V/(a*a)

  Tur(4,1)=(phi*M/4+(T-1)/2)*V_tot
  Tur(4,2)=-U*(phi*M/2+T-1)
  Tur(4,3)=-V*(phi*M/2+T-1)
  Tur(4,4)=phi*M/2+T
 
 !Part 27:
  Do I=1,4
    D(I)=Tur(I,1)*Con(1,J)+Tur(I,2)*Con(2,J)+Tur(I,3)*Con(3,J)+Tur(I,4)*Con(4,J)
  End Do
  Con(1:4,J)=D(1:4) 

 End Do		  
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
