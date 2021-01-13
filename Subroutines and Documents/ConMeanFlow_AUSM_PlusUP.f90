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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_AUSM_PlusUP(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R,F
 Real(8)::Alpha,Beta,Ku,Kp,Const1,GM,Sigma,Temp1,Temp2,Temp3,U,V,Q,Pm,Rho_L,Rho_R,u_L,u_R,&
          v_L,v_R,P_L,P_R,a_L,a_R,Minf,astar_L,astar_R,Utot_L,Utot_R,Ucontra_L,Ucontra_R,&
          a_Face,Ma_L,Ma_R,Ma_Plus_L,Ma_minus_R,Ma_bar,Mo,M_Co,fa,Ma_p,Ma_face,MassFlux,&
          P_plus_L,P_minus_R,P_u,P_face,H_L,H_R,F1,F2,F3,F4,NXX,NYY,DAA,Kapa
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB   
!*********************************************************************************************	
!Part1 :
 Beta   = 0.125
 Ku     = 0.75
 Kp     = 0.25
 Sigma  = 1.0
 Kapa   = 4.0
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

   !Part 11:
    Rho_L = WNP1(1,L)
    U_L   = WNP1(2,L) / Rho_L
    V_L   = WNP1(3,L) / Rho_L    
    P_L   = P(L) 
     
   !Part 12:
    Rho_R = WNP1(1,R)
    U_R   = WNP1(2,R) / Rho_R
    V_R   = WNP1(3,R) / Rho_R 
    P_R   = P(R)
    
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
    Ma_bar =  ( Ucontra_L*Ucontra_L + Ucontra_R*Ucontra_R ) / ( 2*a_face*a_face )
  
   !Part 19:
    M_Co = Kapa*Minf*Minf !max(0.3,0.5*Minf)

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
