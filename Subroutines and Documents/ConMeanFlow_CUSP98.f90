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
!// Developed by: B. Jodeiri, Mechanical Eng., Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_CUSP98(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,P1,P2,L,R,Phi
 Real(8)::U,V,DAA,NXX,NYY,P_L,R_L,U_L,V_L,H_L,P_R,R_R,U_R,V_R,H_R,a_L,a_R,M_L,M_R,F1_R,F2_R,&
          F3_R,F4_R,F1_L,F2_L,F3_L,F4_L,U_bar,V_bar,H_bar,C_bar,Vt_bar,Mm,Part1,Part2,GM,Pm,&
		  Q,Landa_Plus,Landa_Minus,beta,aStrc,Dis_1,Dis_2,Dis_3,Dis_4, F1,F2,F3,F4,Phi_R,Phi_L 
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Integer,Dimension(1:4,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 Phi=1	

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
    NXX = NX(I)
    NYY = NY(I)

   !Part 6:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 7:
    Q  = U*NXX + V*NYY
    Pm = WB(5,I)

   !Part 8:
    F1 = Q*Wb(1,I) 
    F2 = Q*Wb(2,I) + Pm*NXX
    F3 = Q*Wb(3,I) + Pm*NYY
    F4 = Q*Wb(4,I) + Pm*Q

   !Part 9:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
 End Do

!Part 10:
 DO I=NF1+1,NF2

   !Part 11:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 12:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

   !Part 13:
	P_L =  P(L)
    R_L =  WNP1(1,L)
    U_L =  WNP1(2,L)     /R_L
    V_L =  WNP1(3,L)     /R_L
    H_L = (WNP1(4,L)+P_L)/R_L
    
   !Part 14:
	P_R =  P(R)
    R_R =  WNP1(1,R) 
    U_R =  WNP1(2,R)     /R_R
    V_R =  WNP1(3,R)     /R_R
    H_R = (WNP1(4,R)+P_R)/R_R

   !Part 15:
    a_L  = Dsqrt(GM*P_L/R_L)
	a_R  = Dsqrt(GM*P_R/R_R)

	M_L = (U_L*NXX+V_L*NYY) /  a_L
	M_R = (U_R*NXX+V_R*NYY) /  a_R

   !Part 16:
    F1_R = M_R * WNP1(1,R)      * a_R            
    F2_R = M_R * WNP1(2,R)      * a_R + P_R*NXX
    F3_R = M_R * WNP1(3,R)      * a_R + P_R*NYY
    F4_R = M_R *(WNP1(4,R)+P_R) * a_R           

	F1_L = M_L * WNP1(1,L)      * a_L            
    F2_L = M_L * WNP1(2,L)      * a_L + P_L*NXX 
    F3_L = M_L * WNP1(3,L)      * a_L + P_L*NYY
    F4_L = M_L *(WNP1(4,L)+P_L) * a_L            

   !Part 17:
    R_L = DSQRT(R_L)
	R_R = DSQRT(R_R)

    U_bar = ( U_L*R_L+U_R*R_R ) / ( R_L+R_R )
    V_bar = ( V_L*R_L+V_R*R_R ) / ( R_L+R_R )
    H_bar = ( H_L*R_L+H_R*R_R ) / ( R_L+R_R )

    C_bar = DSQRT((GM-1.0d0)*(H_bar-0.5d0*(U_bar*U_bar+V_bar*V_bar)))
   
    Vt_bar = U_bar*NXX+V_bar*NYY
   
   !Part 18:
    Mm=Vt_bar/C_bar
   
   !Part 19:
    Part1 = ((GM+1)*Vt_bar) / (2*GM)
	Part2 = DSQRT(   ( ((GM-1)*Vt_bar)/(2*GM) )**2 + (C_bar*C_bar/GM)   )

    Landa_Plus  = Part1 + Part2
    Landa_Minus = Part1 - Part2

   !Part 20:
    If( Dabs(Mm)>=1.0d0 )Then
     beta=Mm/Dabs(Mm)
    Elseif( Mm>=0.0d0 .And. Mm<1.0d0 ) Then
     beta=Dmax1(0.0d0,(Vt_bar+Landa_Minus)/(Vt_bar-Landa_Minus))
    Else
     beta=-Dmax1(0.0d0,(Vt_bar+Landa_Plus)/(Vt_bar-Landa_Plus))
    End If
   
   !Part 21:
    If( Mm>0.0d0 .And. Mm<1.0d0 .And. beta>0.0d0 ) Then
     aStrc=-(1.0d0+beta)*Landa_Minus
    Elseif( beta==0.0d0 )then
     aStrc=Dabs(Vt_bar)
    Elseif( Mm<0.0d0 .And. Mm>-1.0d0 .And. beta<0.0d0 ) Then
     aStrc=(1.0d0-beta)*Landa_Plus
    Else 
     aStrc=0.0d0
    End If

   !Part 22:
    IF(Phi==1)Then
     Phi_R = WNP1(4,R)+P_R
     Phi_L = WNP1(4,L)+P_L
    ElseIF(Phi==0)Then
     Phi_R = WNP1(4,R)
     Phi_L = WNP1(4,L)
	EndIF

    Dis_1 = 0.5d0* ( aStrc*( WNP1(1,R)-WNP1(1,L) )  +  beta*( F1_R - F1_L ) )
    Dis_2 = 0.5d0* ( aStrc*( WNP1(2,R)-WNP1(2,L) )  +  beta*( F2_R - F2_L ) )
    Dis_3 = 0.5d0* ( aStrc*( WNP1(3,R)-WNP1(3,L) )  +  beta*( F3_R - F3_L ) )
    Dis_4 = 0.5d0* ( aStrc*( Phi_R    -Phi_L     )  +  beta*( F4_R - F4_L ) )

   !Part 23:
    F1 = ( 0.5d0*( F1_R+F1_L ) - Dis_1 ) * DAA 
    F2 = ( 0.5d0*( F2_R+F2_L ) - Dis_2 ) * DAA 
    F3 = ( 0.5d0*( F3_R+F3_L ) - Dis_3 ) * DAA 
    F4 = ( 0.5d0*( F4_R+F4_L ) - Dis_4 ) * DAA 

   !Part 24:
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4

   !Part 25:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4

 End Do
!*********************************************************************************************
 End
!###########################################################################################

