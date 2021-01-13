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
 Subroutine Limiter3D(Dim,NC,NF1,NF2,IDS,GM,XC,YC,ZC,FaceType,WNP1,P,WB,GWNP1,Limit,X,Y,Z,NF)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,IDS,GM,XC,YC,ZC,FaceType,WNP1,P,WB,GWNP1,X,Y,Z,NF
 Intent(Out  )::Limit

 Integer::Dim,I,J,K,NC,NF1,NF2,ME,NE,NFacePnt
 Real(8)::Del_Minus,Del_Plus,Phi_ik,U_ik,GM,Min_U_min,Max_U_max,eps,Temp1,Temp2,Temp3,Temp4,&
          Num,Denom,C,Utot,DelU,DX,DY,DZ,Fc,Phi
 Integer,Dimension(1:6,1:Dim)::IDS 
 Real(8),Dimension(1:Dim)::Ma,MaxMa,U_max,U_min,XC,YC,ZC,U,P
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5,1:Dim)::WNP1,Limit
 Real(8),Dimension(1:3,1:5,1:Dim)::GWNP1
 
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:4)::PG
 Real(8)::XM_Face,YM_Face,zM_Face
 
 Real(8),Dimension(1:Dim)::X,Y,Z
 
 Real(8)::RRR,R_G,U_G,V_G,P_G,M_G
 Integer::NF !Number of Faces Constructing Mesh
 
 Real(8)::EPS_BRT
!*********************************************************************************************
!Part 1:   
 Limit = 1.0  

 Do J=1,5  

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
        U(I) = WNP1(4,I)/WNP1(1,I)                               
     End Do
    ElseIF(J==5)Then
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
        
           RRR = WB(4,I)/WB(1,I)
           
       ElseIF(J==5)Then
        
           RRR = WB(6,I)
           
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
       
       !Part 17:
	   DX = XM_Face - XC(ME)
       DY = YM_Face - YC(ME)
       DZ = ZM_Face - ZC(ME)
      
      !Part 18:
       Fc = U(ME)
       Phi = 1.0
       
       Call Recons2Ord3D(Dim,J,ME,Fc,GWNP1,Phi,DX,DY,DZ,U_ik)
     
      !Part 19:
       Del_Minus = (U_ik-U(ME))
       
       IF (Del_Minus > 0) Del_Plus = U_Max(ME) - U(ME)
       IF (Del_Minus < 0) Del_Plus = U_Min(ME) - U(ME)
                   
      !Part 20:
       Num = ( ( Del_Plus**2.0 + eps )* Del_Minus ) + ( 2.0*( Del_Minus**2.0 )* Del_Plus )
       Denom = Del_Minus * ( Del_Plus**2.0 +  2.0*( Del_Minus**2.0 ) + Del_Minus*Del_Plus + eps )


       IF (Del_Minus == 0.0) THEN
       Phi_ik=1.0
       ELSE
       Phi_ik =   Num / Denom 
       END IF 
                 
      !Part 21:
       IF ( Phi_ik < Limit(J,ME) ) Limit(J,ME)= Phi_ik 
             
		    
      !Part 22:
	   DX = XM_Face - XC(ME)
       DY = YM_Face - YC(ME)
       DZ = ZM_Face - ZC(ME)
       
       Fc = U(NE)
       Phi = 1.
       Call Recons2Ord3D(Dim,J,NE,Fc,GWNP1,Phi,DX,DY,DZ,U_ik) 
     
       
       Del_Minus = (U_ik-U(NE))

       IF (Del_Minus > 0) Del_Plus = U_Max(NE) - U(NE)
       IF (Del_Minus < 0) Del_Plus = U_Min(NE) - U(NE)

       
       Num = ( ( Del_Plus**2.0 + eps )* Del_Minus ) + ( 2.0*( Del_Minus**2.0 )* Del_Plus )
       Denom = Del_Minus * ( Del_Plus**2.0 +  2.0*( Del_Minus**2.0 ) + Del_Minus*Del_Plus + eps )

       
       IF (Del_Minus == 0.0) THEN
       Phi_ik=1.0
       ELSE
       Phi_ik =   Num / Denom 
       END IF 
              
       
       IF ( Phi_ik < Limit(J,NE) ) Limit(J,NE)= Phi_ik
          
    End Do

    !Part 23:
    Do I=NF2+1, NF 
       
      !Part 24:
       ME = IDS(1,I)
              
       !Part 25:
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
       
       !Part 26:
	   DX = XM_Face - XC(ME)
       DY = YM_Face - YC(ME)
       DZ = ZM_Face - ZC(ME)
      
      !Part 27:
       Fc = U(ME)
       Phi = 1.
       Call Recons2Ord3D(Dim,J,ME,Fc,GWNP1,Phi,DX,DY,DZ,U_ik)
     
      !Part 28:
       Del_Minus = (U_ik-U(ME))
       
       IF (Del_Minus > 0) Del_Plus = U_Max(ME) - U(ME)
       IF (Del_Minus < 0) Del_Plus = U_Min(ME) - U(ME)
                   
      !Part 29:
       Num = ( ( Del_Plus**2.0 + eps )* Del_Minus ) + ( 2.0*( Del_Minus**2.0 )* Del_Plus )
       Denom = Del_Minus * ( Del_Plus**2.0 +  2.0*( Del_Minus**2.0 ) + Del_Minus*Del_Plus + eps )

       IF (Del_Minus == 0.0) THEN
       Phi_ik=1.0
       ELSE
       Phi_ik =   Num / Denom 
       END IF 
                   
      !Part 30:
       IF ( Phi_ik < Limit(J,ME) ) Limit(J,ME)= Phi_ik 
             
    End Do
     
 End Do    

!*********************************************************************************************
 End
!##############################################################################################
 
 
