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
 Subroutine Limiter(Dim,NC,NF1,NF2,IDS,GM,XC,YC,WNP1,P,WB,GWNP1,Limit,X,Y,NF)
 Implicit None !V2 az name file bardashte shod
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,IDS,GM,XC,YC,WNP1,P,WB,GWNP1,X,Y,NF
 Intent(Out  )::Limit

 Integer::Dim,I,J,NC,NF1,NF2,ME,NE
 Real(8)::Del_Minus,Del_Plus,Phi_ik,U_ik,GM,Min_U_min,Max_U_max,eps,Temp1,Temp2,Temp3,Temp4,&
          Num,Denom,C,Utot,DelU,DX,DY,Fc,Phi
 Integer,Dimension(1:4,1:Dim)::IDS 
 Real(8),Dimension(1:Dim)::Ma,MaxMa,U_max,U_min,XC,YC,U,P
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:4,1:Dim)::WNP1,Limit
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 
 INTEGER::P1,P2
 Real(8)::X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG
 Real(8),Dimension(1:Dim)::X,Y
 
 Real(8)::RRR,R_G,U_G,V_G,P_G,M_G
 Integer::NF !Number of Faces Constructing Mesh
 
 Real(8)::EPS_BRT
!*********************************************************************************************
!Part 1:   
 Limit = 1.0  

 Do J=1, 4  

   !Part 2: 
    IF(J==1)Then
     Do I=1,NC 
        U(I) = WNP1(1,I)                                
     End Do
	ElseIF(J==2)Then
     Do I=1,NC 
        U(I) = WNP1(2,I)/WNP1(1,I)                                
     End Do
	ElseIF(J==3)Then
     Do I=1,NC 
        U(I) = WNP1(3,I)/WNP1(1,I)                               
     End Do
	ElseIF(J==4)Then
     Do I=1,NC 
        U(I) = P(I)                               
     End Do
    EndIF

    !Part 3:
    Do I=1, NC 
       U_min(I) = U(I)
       U_max(I) = U(I)                 
    End Do
    
   !Part 4:
    Do I=NF1+1, NF2 
    
      !Part 5:
       ME = IDS(1,I)
       NE = IDS(2,I)
       
      !Part 6:
       Temp1 = U(ME)
       Temp2 = U(NE)
       
       IF (Temp2 > U_max(ME) ) U_max(ME) = Temp2
       IF (Temp2 < U_min(ME) ) U_min(ME) = Temp2
       IF (Temp1 > U_max(NE) ) U_max(NE) = Temp1
       IF (Temp1 < U_min(NE) ) U_min(NE) = Temp1
           
    End Do
    
   
    !Part 7:
    Do I=NF2+1, NF 
    
      !Part 8:
       ME = IDS(1,I)
        
       !Part 9:
       IF(J==1)Then
        
           RRR = WB(1,I)
           
       ElseIF(J==2)Then
        
           RRR = WB(2,I)/WB(1,I)
           
       ElseIF(J==3)Then
        
           RRR = WB(3,I)/WB(1,I)
           
       ElseIF(J==4)Then
        
           RRR = WB(5,I)
           
       EndIF
       
       !Part 10:
       R_G = 2.0*RRR-U(ME)
       
       !Part 11:
       Temp1 = U(ME)
       Temp2 = R_G
    
       IF (Temp2 > U_max(ME) ) U_max(ME) = Temp2
       IF (Temp2 < U_min(ME) ) U_min(ME) = Temp2
       
          
    End Do
    
    

   !Part 12:
    Max_U_max =-1000.0
    Min_U_min = 1000.0
    Do I=1, NC     
       IF (U_max(I) > Max_U_max)  Max_U_max = U_max(I)
       IF (U_min(I) < Min_U_min)  Min_U_min = U_min(I)
    End Do
    
    !Part 13:
    eps = 0.05* ( Max_U_max - Min_U_min )
    eps = eps*eps

   !Part 14:
    Do I=NF1+1,NF2
       
      !Part 15:
       ME = IDS(1,I)
       NE = IDS(2,I)   
       
       !Part 16:
       P1=IDS(3,I)          
       P2=IDS(4,I)          
                            
       X_P1=X(P1)
       Y_P1=Y(P1)
       
       X_P2=X(P2)
       Y_P2=Y(P2)
       
       XM_EDG=0.5*(X_P1+X_P2)
       YM_EDG=0.5*(Y_P1+Y_P2)
       
       !Part 17:
	   DX = XM_EDG - XC(ME)
       DY = YM_EDG - YC(ME)
      
      !Part 18:
       Fc = U(ME)
       Phi = 1.
       Call Recons2Ord(Dim,J,ME,Fc,GWNP1,Phi,DX,DY,U_ik)
     
      !Part 19:
       Del_Minus = (U_ik-U(ME))
       
       IF (Del_Minus > 0) Del_Plus = U_Max(ME) - U(ME)
       IF (Del_Minus < 0) Del_Plus = U_Min(ME) - U(ME)
                   
      !Part 20:
       Num = ( ( Del_Plus**2.0 + eps )* Del_Minus ) + ( 2.0*( Del_Minus**2.0 )* Del_Plus )
       Denom = Del_Minus * ( Del_Plus**2.0 +  2.0*( Del_Minus**2.0 ) + Del_Minus*Del_Plus + eps )
       Phi_ik =   Num / Denom 
       
       IF (Del_Minus == 0.0) Phi_ik=1.0
                 
      !Part 21:
       IF ( Phi_ik < Limit(J,ME) ) Limit(J,ME)= Phi_ik 
             
		    
      !Part 22:
	   DX = XM_EDG - XC(NE)
       DY = YM_EDG - YC(NE)
		
       
       Fc = U(NE)
       Phi = 1.
       Call Recons2Ord(Dim,J,NE,Fc,GWNP1,Phi,DX,DY,U_ik) 
     
       
       Del_Minus = (U_ik-U(NE))

       IF (Del_Minus > 0) Del_Plus = U_Max(NE) - U(NE)
       IF (Del_Minus < 0) Del_Plus = U_Min(NE) - U(NE)

       
       Num = ( ( Del_Plus**2.0 + eps )* Del_Minus ) + ( 2.0*( Del_Minus**2.0 )* Del_Plus )
       Denom = Del_Minus * ( Del_Plus**2.0 +  2.0*( Del_Minus**2.0 ) + Del_Minus*Del_Plus + eps )
       Phi_ik =   Num / Denom
       
       IF (Del_Minus == 0.0) Phi_ik=1.0
              
       
       IF ( Phi_ik < Limit(J,NE) ) Limit(J,NE)= Phi_ik
          
    End Do

    !Part 23:
    Do I=NF2+1, NF 
       
      !Part 24:
       ME = IDS(1,I)
              
       !Part 25:
       P1=IDS(3,I)          
       P2=IDS(4,I)          
                            
       X_P1=X(P1)
       Y_P1=Y(P1)
       
       X_P2=X(P2)
       Y_P2=Y(P2)
       
       XM_EDG=0.5*(X_P1+X_P2)
       YM_EDG=0.5*(Y_P1+Y_P2)
       
       !Part 26:
	   DX = XM_EDG - XC(ME)
       DY = YM_EDG - YC(ME)
      
      !Part 27:
       Fc = U(ME)
       Phi = 1.
       Call Recons2Ord(Dim,J,ME,Fc,GWNP1,Phi,DX,DY,U_ik)
     
      !Part 28:
       Del_Minus = (U_ik-U(ME))
       
       IF (Del_Minus > 0) Del_Plus = U_Max(ME) - U(ME)
       IF (Del_Minus < 0) Del_Plus = U_Min(ME) - U(ME)
                   
      !Part 29:
       Num = ( ( Del_Plus**2.0 + eps )* Del_Minus ) + ( 2.0*( Del_Minus**2.0 )* Del_Plus )
       Denom = Del_Minus * ( Del_Plus**2.0 +  2.0*( Del_Minus**2.0 ) + Del_Minus*Del_Plus + eps )
       Phi_ik =   Num / Denom 
       
       IF (Del_Minus == 0.0) Phi_ik=1.0
                   
      !Part 30:
       IF ( Phi_ik < Limit(J,ME) ) Limit(J,ME)= Phi_ik 
             
    End Do
     
 End Do   
 

!*********************************************************************************************
 End
!##############################################################################################
 
 
