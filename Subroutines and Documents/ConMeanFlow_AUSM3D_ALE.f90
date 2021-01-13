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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_AUSM3D_ALE(Dim,NC,NF1,NF2,NF,GM,IDS,Nx,Ny,Nz,DA,WNP1,WB,P,GF,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,GM,IDS,Nx,Ny,Nz,DA,WNP1,WB,P,GF
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,P1,P2,L,R
 Real(8)::U,V,W,DX,DY,F1,F2,F3,F4,F5,Q,Ro,RU,RV,RW,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,&
          P_Minus,Mm,Pm,DL,Nxx,Nyy,Nzz,DS
 Real(8),Dimension(1:5,1:Dim)::WNP1,Con
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::Nx,Ny,Nz,DA,P,GF
 Integer,Dimension(1:6,1:Dim)::IDS
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
    Con(5,I) = 0.0
 End Do
 
!Part 2:
 DO I=NF2+1,NF
     
   !Part 3:
    ME = IDS(1,I)
    
   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)
    
	NXX = NX(I)    
	NYY = NY(I)    
	NZZ = NZ(I)

   !Part 5:
    Q  = U*NXX+V*NYY+W*NZZ
    Pm = WB(6,I)
 
   !Part 6:
    F1 = ( Q - GF(I) ) * WB(1,I) 
    F2 = ( Q - GF(I) ) * WB(2,I) + Pm*NXX
    F3 = ( Q - GF(I) ) * WB(3,I) + Pm*NYY
    F4 = ( Q - GF(I) ) * WB(4,I) + Pm*NZZ
    F5 = ( Q - GF(I) ) * WB(5,I) + Pm*Q
    
   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
    Con(5,ME) = Con(5,ME) + F5
 End Do
 
!Part 8:
 DO I=NF1+1,NF2
     
   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DS  = DA(I)
	NXX = NX(I) / DS  
	NYY = NY(I) / DS     
	NZZ = NZ(I) / DS 

   !Part 11:
    a_L  = Dsqrt(DABS(GM*P(L)/Wnp1(1,L)))
	a_R  = Dsqrt(DABS(GM*P(R)/Wnp1(1,R)))

	M_L = ( (Wnp1(2,L)*NXX + Wnp1(3,L)*NYY + Wnp1(4,L)*NZZ)/Wnp1(1,L) - (GF(I)/DS) ) / a_L
    M_R = ( (Wnp1(2,R)*NXX + Wnp1(3,R)*NYY + Wnp1(4,R)*NZZ)/Wnp1(1,R) - (GF(I)/DS) ) / a_R
    
   !Part 12:
    IF(Dabs(M_L)>1.)Then
       M_Plus = 0.5*(M_L+Dabs(M_L))
       P_Plus = 0.5*(M_L+Dabs(M_L))/(M_L)   !!!!???? /Dabs(M_L)   
    Else
       M_Plus = 0.25*(M_L+1.)*(M_L+1.)
       P_Plus = 0.25*(M_L+1.)*(M_L+1.)*(2.-M_L)
    End If
    
   !Part 13:
    IF(Dabs(M_R)>1.)Then
        M_Minus = 0.5*(M_R-Dabs(M_R))
        P_Minus = 0.5*(M_R-Dabs(M_R))/(M_R)    !!!!????  /Dabs(M_R)
    Else
        M_Minus =-0.25*(M_R-1.)*(M_R-1.)
        P_Minus = 0.25*(M_R-1.)*(M_R-1.)*(2.+M_R)
    End If
    
   !Part 14:
    Mm = M_Plus+M_Minus
    Pm = P(L)*P_Plus + P(R)*P_Minus
    
   !Part 15:
    If(Mm<=0.)Then
        Ro = Wnp1(1,R)       * a_R
        RU = Wnp1(2,R)       * a_R
        RV = Wnp1(3,R)       * a_R
        RW = Wnp1(4,R)       * a_R
        RH =(Wnp1(5,R)+P(R)) * a_R
    Else
        Ro = Wnp1(1,L)       * a_L
        RU = Wnp1(2,L)       * a_L
        RV = Wnp1(3,L)       * a_L
        RW = Wnp1(4,L)       * a_L
        RH =(Wnp1(5,L)+P(L)) * a_L
    End If
    
   !Part 16:
    F1 = ( Mm * Ro          ) * DS
    F2 = ( Mm * RU + Pm*NXX ) * DS
    F3 = ( Mm * RV + Pm*NYY ) * DS
    F4 = ( Mm * RW + Pm*NZZ ) * DS
    F5 = ((Mm * RH          ) * DS)+(Pm * GF(I))
    
   !Part 17:
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4
    Con(5,L) = Con(5,L) + F5
    
   !Part 18:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    Con(5,R) = Con(5,R) - F5
    
 End Do
!*********************************************************************************************
 End
!###########################################################################################
