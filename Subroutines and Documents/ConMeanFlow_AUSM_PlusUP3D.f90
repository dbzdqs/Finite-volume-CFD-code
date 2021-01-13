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
 Subroutine ConMeanFlow_AUSM_PlusUP3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Minf,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Minf,WNP1,WB,P
 Intent(Out  )::Con

 Integer::I,ME,NE,L,R,F
 Real(8)::Alpha,Beta,Ku,Kp,Const1,GM,Sigma,Temp1,Temp2,Temp3,U,V,W,Q,Pm,Rho_L,Rho_R,u_L,u_R,&
          v_L,v_R,w_L,w_R,P_L,P_R,a_L,a_R,Minf,astar_L,astar_R,Utot_L,Utot_R,Ucontra_L,Ucontra_R,&
          a_Face,Ma_L,Ma_R,Ma_Plus_L,Ma_minus_R,Ma_bar,Mo,M_Co,fa,Ma_p,Ma_face,MassFlux,&
          P_plus_L,P_minus_R,P_u,P_face,H_L,H_R,F1,F2,F3,F4,F5,NXX,NYY,NZZ,DAA,Kapa
!--------------------------------------------------------------------------------------------
 Integer::Dim !Maximum Dimension of Arrays
 Integer::NC !Number of Existing Cells
 Integer::NF1,NF2 !Index of 1st and last Non-Boundary Faces
 Integer::NF !Number of Faces Constructing Mesh !Index of Last Face of Mesh
 Integer,Dimension(1:6,1:Dim)::IDS !Information of Grid Data Structure  
 Real(8),Dimension(1:Dim)::NX,NY,NZ !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA !Area of each Face
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:5,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(1:6,1:Dim)::WB !Conservative Values and Pressure at Boundary Faces
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
    Con(5,I) = 0.0
 End Do  
   
!Part 3:
 DO I=NF2+1,NF

   !Part 4:
    ME = IDS(1,I)

   !Part 5:
    U = WB(2,I) / WB(1,I)
    V = WB(3,I) / WB(1,I)
    W = WB(4,I) / WB(1,I)

   !Part 6:
	Nxx = Nx(I)    
	Nyy = Ny(I)    
	Nzz = Nz(I)
	
    Q  = U*Nxx+V*Nyy+W*Nzz
    Pm = WB(6,I)

   !Part 7:
    F1 = Q * WB(1,I) 
    F2 = Q * WB(2,I)     + Pm*Nxx
    F3 = Q * WB(3,I)     + Pm*Nyy
    F4 = Q * WB(4,I)     + Pm*Nzz
    F5 = Q *(WB(5,I)+Pm)

   !Part 8:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
    Con(5,ME) = Con(5,ME) + F5  

 End Do

!Part 9:
 DO I=NF1+1,NF2

   !Part 10:
    L = IDS(1,I)
    R = IDS(2,I)

    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA
    NZZ = NZ(I)/DAA

   !Part 11:
    Rho_L = WNP1(1,L)
    U_L   = WNP1(2,L) / Rho_L
    V_L   = WNP1(3,L) / Rho_L  
    W_L   = WNP1(4,L) / Rho_L  
    P_L   = P(L) 
     
   !Part 12:
    Rho_R = WNP1(1,R)
    U_R   = WNP1(2,R) / Rho_R
    V_R   = WNP1(3,R) / Rho_R 
    W_R   = WNP1(4,R) / Rho_R 
    P_R   = P(R)
    
   !Part 13:
   !Left
    Utot_L  = u_L*u_L + v_L*v_L + w_L*w_L
    H_L     = ( (GM * P_L / Rho_L) / (GM - 1.0) ) + (0.5 * Utot_L)
    astar_L = Dsqrt(H_L  * Const1)

   !Right
    Utot_R  = u_R*u_R + v_R*v_R + w_R*w_R 
    H_R     = ( ((GM * P_R / Rho_R) / (GM - 1.0) ) + (0.5 * Utot_R))

    astar_R = Dsqrt(H_R  * Const1)    

   !Part 14:
    Ucontra_L = u_L*NXX + v_L*NYY  + w_L*NZZ 
    Ucontra_R = u_R*NXX + v_R*NYY  + w_R*NZZ
    
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
     F4 = ( MassFlux * w_L + P_Face * NZZ ) * DAA
     F5 = ( MassFlux * H_L                ) * DAA
    Else
     F1 = ( MassFlux                      ) * DAA 
     F2 = ( MassFlux * u_R + P_Face * NXX ) * DAA
     F3 = ( MassFlux * v_R + P_Face * NYY ) * DAA
     F4 = ( MassFlux * w_R + P_Face * NZZ ) * DAA
     F5 = ( MassFlux * H_R                ) * DAA
    End IF          

   !Part 27:
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4
    Con(5,L) = Con(5,L) + F5
	    
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    Con(5,R) = Con(5,R) - F5      

 End Do 
!********************************************************************************************* 
 End
!###########################################################################################
