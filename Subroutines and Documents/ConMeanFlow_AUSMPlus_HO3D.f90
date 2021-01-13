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
!// Date: Dec., 05, 2016                                                                   //!
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
 Subroutine ConMeanFlow_AUSM_Plus_HO3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,FaceType,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Zc,Limit,Con,X,Y,Z)
 Implicit none
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,FaceType,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Zc,X,Y,Z
 Intent(Out  )::Con

 Integer::Dim,I,K,NC,NF1,NF2,NF,ME,NE,R,L,NFacePnt
 Real(8)::Alpha,Beta,Const1,GM,U,V,W,Q,Pm,Rho_L,Rho_R,U_L,U_R,V_L,V_R,W_L,W_R,P_L,P_R,a_L,a_R,astar_L,&
          astar_R,Utot_L,Utot_R,Ucontra_L,Ucontra_R,a_Face,Ma_L,Ma_R,Ma_Plus_L,Ma_minus_R,Phi,&
          Ma_face,MassFlux,P_plus_L,P_minus_R,P_face,H_L,H_R,F1,F2,F3,F4,F5,DAA,NXX,NYY,NZZ,Ro,Pr,DX,DY,DZ
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P,Xc,Yc,Zc
 Real(8),Dimension(1:5,1:Dim)::Con,WNP1,Limit
 Real(8),Dimension(1:6,1:Dim)::WB    
 Real(8),Dimension(1:3,1:5,1:Dim)::GWNP1
 
  Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:4)::PG
 Real(8)::XM_Face,YM_Face,zM_Face
 
 Real(8),Dimension(1:Dim)::X,Y,Z
 INTEGER::P1,P2,P3,P4
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,X_P3,Y_P3,X_P4,Y_P4,XM_EDG,YM_EDG
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
    Q = U*NX(I) + V*NY(I) + W*NZ(I)        
    Pm = WB(6,I)

   !Part 7:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I) 
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q * Wb(4,I) + Pm*NZ(I)
    F5 = Q *(Wb(5,I)+Pm)

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
      XM_Face=0.0
      YM_Face=0.0
      ZM_Face=0.0  
       
     NFacePnt=FaceType(I)
 
   DO K=1,NFacePnt   
       PG(K) = IDS(K+2,I)
       XM_Face=XM_Face+X(PG(K))
       YM_Face=YM_Face+Y(PG(K))
       ZM_Face=ZM_Face+Z(PG(K))
   END DO
    
       XM_Face = XM_Face / NFacePnt
       YM_Face = YM_Face / NFacePnt
       ZM_Face = ZM_Face / NFacePnt
       

      !Part 12:
       Ro = WNP1(1,L)
       U  = WNP1(2,L) / Ro 
       V  = WNP1(3,L) / Ro  
       W  = WNP1(4,L) / Ro           
       Pr = P(L) 
    
       DX = XM_Face - XC(L)
       DY = YM_Face - YC(L)     
       DZ = ZM_Face - ZC(L)     
       
       Phi=Limit(1,L)
	   Call Recons2Ord3D(Dim,1,L,Ro,GWNP1,Phi,DX,DY,DZ,Ro) 
       Phi=Limit(2,L)        
	   Call Recons2Ord3D(Dim,2,L,U,GWNP1,Phi,DX,DY,DZ,U)
       Phi=Limit(3,L)         
	   Call Recons2Ord3D(Dim,3,L,V,GWNP1,Phi,DX,DY,DZ,V)
       Phi=Limit(4,L)
       Call Recons2Ord3D(Dim,4,L,W,GWNP1,Phi,DX,DY,DZ,W)
       Phi=Limit(5,L)         
	   Call Recons2Ord3D(Dim,5,L,Pr,GWNP1,Phi,DX,DY,DZ,Pr)

       Rho_L = Ro
       U_L   = U 
       V_L   = V  
       W_L   = W 
       P_L   = Pr

      !Part 13:
       Ro = WNP1(1,R)
       U  = WNP1(2,R) / Ro 
       V  = WNP1(3,R) / Ro  
       W  = WNP1(4,R) / Ro           
       Pr = P(R) 
    
       DX = XM_Face - XC(R)
       DY = YM_Face - YC(R)     
       DZ = ZM_Face - ZC(R)  

       Phi=Limit(1,R)
	   Call Recons2Ord3D(Dim,1,R,Ro,GWNP1,Phi,DX,DY,DZ,Ro)
       Phi=Limit(2,R)         
	   Call Recons2Ord3D(Dim,2,R,U,GWNP1,Phi,DX,DY,DZ,U)
       Phi=Limit(3,R)         
	   Call Recons2Ord3D(Dim,3,R,V,GWNP1,Phi,DX,DY,DZ,V)
       Phi=Limit(4,R) 
       Call Recons2Ord3D(Dim,4,R,W,GWNP1,Phi,DX,DY,DZ,W) 
       Phi=Limit(5,R)        
	   Call Recons2Ord3D(Dim,5,R,Pr,GWNP1,Phi,DX,DY,DZ,Pr)
  
       Rho_R = Ro
       U_R   = U 
       V_R   = V 
       W_R   = W 
       P_R   = Pr

      !Part 14:
       Utot_L  = U_L*U_L + V_L*V_L + W_L*W_L
       H_L     = ( (  GM * P_L / Rho_L  ) / ( GM - 1.0 ) ) + ( 0.5 * Utot_L )
       astar_L = Dsqrt(H_L  * Const1)

       Utot_R  =  U_R*U_R + V_R*V_R + W_R*W_R
       H_R     = ( ( GM * P_R / Rho_R ) / ( GM - 1.0 ) ) + ( 0.5 * Utot_R ) 
       astar_R = Dsqrt(H_R  * Const1)   
   
      !Part 15:
       Ucontra_L = U_L*NXX + V_L*NYY + W_L*NZZ
       Ucontra_R = U_R*NXX + V_R*NYY + W_R*NZZ

      !Part 16:
       a_L = ( astar_L*astar_L ) / DMax1( astar_L , Ucontra_L )
       a_R = ( astar_R*astar_R ) / DMax1( astar_R ,-Ucontra_R )

       a_face = DMax1( a_L , a_R )

      !Part 17:    
       Ma_L = Ucontra_L / a_face
       Ma_R = Ucontra_R / a_face    
      
      !Part 18: 
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

      !Part 19:  
       Ma_face = Ma_plus_L + Ma_minus_R       
  
      !Part 20: 
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

      !Part 21:
       P_Face = P_plus_L * P_L + P_minus_R * P_R

      !Part 22:    
       IF ( Ma_face >= 0.0 )Then 
        MassFlux = a_face * Ma_face * Rho_L
       Else
        MassFlux = a_face * Ma_face * Rho_R
       End IF  
    
      !Part 23:
       IF ( MassFlux >= 0.0 )Then

 
        F1 = ( MassFlux                      ) * DAA
        F2 = ( MassFlux * U_L + P_Face * NXX ) * DAA
        F3 = ( MassFlux * V_L + P_Face * NYY ) * DAA
        F4 = ( MassFlux * W_L + P_Face * NZZ ) * DAA       
        F5 = ( MassFlux * H_L                ) * DAA
       Else
        F1 = ( MassFlux                      ) * DAA
        F2 = ( MassFlux * U_R + P_Face * NXX ) * DAA
        F3 = ( MassFlux * V_R + P_Face * NYY ) * DAA
        F4 = ( MassFlux * W_R + P_Face * NZZ ) * DAA
        F5 = ( MassFlux * H_R                ) * DAA
       End IF   
   
      !Part 24:
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
       
     if( R==1)Then
 print*,'NE', Con(1,R) 
 pause
 endif
      

 End Do
 
 !print*,'=========='

!*********************************************************************************************
 End
!###########################################################################################
