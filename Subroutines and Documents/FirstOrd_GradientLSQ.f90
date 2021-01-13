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
 Subroutine FirstOrd_GradientLSQ(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1,XC,YC,X,Y,R_LSQ)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,XC,YC,X,Y,R_LSQ
 Intent(Out  )::GWNP1

 Integer::Dim,I,J,NC,NF1,NF2,NF,ME,NE
 Real(8)::NXX,NYY,R,U,V,Pr,R_NX,U_NX,V_NX,P_NX,R_NY,U_NY,V_NY,P_NY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,A  ,P    
 Real(8),Dimension(1:4,1:Dim)::WNP1      
 Real(8),Dimension(1:5,1:Dim)::WB    
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 
 Real(8),Dimension(1:Dim)::X,Y,XC,YC
 Real(8),Dimension(1:3,1:Dim)::R_LSQ
 Real(8)::X_L,Y_L,X_R,Y_R,X_P1,X_P2,Y_P1,Y_P2,XM_EDG,YM_EDG,DELX_L,DELY_L,DELX_R,DELY_R
 
 Real(8)::R11_L,R12_L,R22_L,R11_R,R12_R,R22_R,WX_L,WX_R,WY_L,WY_R
 Real(8)::R_ME,U_ME,V_ME,P_ME,R_G,U_G,V_G,P_G
 Real(8)::R_L,U_L,V_L,P_L,R_R,U_R,V_R,P_R
 Integer::P1,P2

 Real(8)::XP_L,YP_L,SLOPE,CNST

!*********************************************************************************************	
!Part 1:    
 GWNP1=0.0
 

!Part 2:    
DO I=NF2+1,NF   

   !Part 3: 
    ME = IDS(1,I) 
   
	NXX = NX(I)
    NYY = NY(I)

   !Part 4:
	R = WB(1,I)
    U = WB(2,I) / R
    V = WB(3,I) / R 
    Pr= WB(5,I)
    
    !Part 5:
    R_ME = WNP1(1,ME)
    U_ME = WNP1(2,ME)/WNP1(1,ME)
    V_ME = WNP1(3,ME)/WNP1(1,ME)
    P_ME = P(ME)
    
    !Part 6:
    X_L=XC(ME)
    Y_L=YC(ME)
    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    
    !Part 7:
    IF ((X_P2-X_P1)/=0.0 .AND. (Y_P2-Y_P1)/=0) THEN
        
        SLOPE = (Y_P2-Y_P1)/(X_P2-X_P1)
        CNST = Y_P1-SLOPE*X_P1
        
        XP_L = (2.0*(Y_L-CNST)+X_L*((1.0/SLOPE)-SLOPE))/(SLOPE+(1.0/SLOPE))
        YP_L = Y_L-(1.0/SLOPE)*(XP_L-X_L)
        
    !Part 8:    
    ELSE IF ((X_P2-X_P1)==0.0) THEN
        
        XP_L = 2.0*X_P1-X_L
        YP_L = Y_L
    
    !Part 9:
    ELSE IF ((Y_P2-Y_P1)==0.0) THEN
        
        XP_L = X_L
        YP_L = 2.0*Y_P1-Y_L
        
    END IF
    
    !Part 10:
    DELX_L = (XP_L-X_L)
    DELY_L = (YP_L-Y_L)
    
    !Part 11:
    R11_L = R_LSQ(1,ME)
    R12_L = R_LSQ(2,ME)
    R22_L = R_LSQ(3,ME)
    
    !Part 12:
    WX_L = DELX_L/(R11_L**2)-(R12_L/(R11_L*(R22_L**2)))*(DELY_L-DELX_L*(R12_L/R11_L))
    WY_L = (1.0/(R22_L**2))*(DELY_L-DELX_L*(R12_L/R11_L))
    
    !!GHOST PRIMITIVE VARIABLES
    !Part 13:
    R_G = 2.0*R-R_ME
    U_G = 2.0*U-U_ME
    V_G = 2.0*V-V_ME
    P_G = 2.0*Pr-P_ME
    
    !Part 14:
    GWNP1(1,1,ME)=GWNP1(1,1,ME) + WX_L*(R_G-R_ME)
    GWNP1(2,1,ME)=GWNP1(2,1,ME) + WY_L*(R_G-R_ME)
    
    GWNP1(1,2,ME)=GWNP1(1,2,ME) + WX_L*(U_G-U_ME)
    GWNP1(2,2,ME)=GWNP1(2,2,ME) + WY_L*(U_G-U_ME) 
      
    GWNP1(1,3,ME)=GWNP1(1,3,ME) + WX_L*(V_G-V_ME)
    GWNP1(2,3,ME)=GWNP1(2,3,ME) + WY_L*(V_G-V_ME)
    
    GWNP1(1,4,ME)=GWNP1(1,4,ME) + WX_L*(P_G-P_ME) 
    GWNP1(2,4,ME)=GWNP1(2,4,ME) + WY_L*(P_G-P_ME)
     
END DO

!Part 15:
DO I=NF1+1,NF2
     
    !Part 16:
	 ME = IDS(1,I)
     NE = IDS(2,I)
     
     !Part 17:
     R_L = WNP1(1,ME)
     U_L = WNP1(2,ME)/WNP1(1,ME)
     V_L = WNP1(3,ME)/WNP1(1,ME)
     P_L = P(ME)
     
     R_R = WNP1(1,NE)
     U_R = WNP1(2,NE)/WNP1(1,NE)
     V_R = WNP1(3,NE)/WNP1(1,NE)
     P_R = P(NE)
     
     !Part 18:
     X_L=XC(ME)
     Y_L=YC(ME)
     
     X_R=XC(NE)
     Y_R=YC(NE)
     
     !Part 19:
     DELX_L = X_R-X_L
     DELY_L = Y_R-Y_L
         
     DELX_R = X_L-X_R
     DELY_R = Y_L-Y_R
       
     !Part 20:
     R11_L = R_LSQ(1,ME)
     R12_L = R_LSQ(2,ME)
     R22_L = R_LSQ(3,ME)
     
     !Part 21:
     WX_L = DELX_L/(R11_L**2)-(R12_L/(R11_L*(R22_L**2)))*(DELY_L-DELX_L*(R12_L/R11_L))
     WY_L = (1.0/(R22_L**2))*(DELY_L-DELX_L*(R12_L/R11_L))
      

      
     !Part 22:
     GWNP1(1,1,ME)=GWNP1(1,1,ME) + WX_L*(R_R-R_L) 
     GWNP1(2,1,ME)=GWNP1(2,1,ME) + WY_L*(R_R-R_L)            
              
     GWNP1(1,2,ME)=GWNP1(1,2,ME) + WX_L*(U_R-U_L)
     GWNP1(2,2,ME)=GWNP1(2,2,ME) + WY_L*(U_R-U_L)             
               
     GWNP1(1,3,ME)=GWNP1(1,3,ME) + WX_L*(V_R-V_L)
     GWNP1(2,3,ME)=GWNP1(2,3,ME) + WY_L*(V_R-V_L)             
              
     GWNP1(1,4,ME)=GWNP1(1,4,ME) + WX_L*(P_R-P_L)
     GWNP1(2,4,ME)=GWNP1(2,4,ME) + WY_L*(P_R-P_L)            
     
     !Part 23:
     R11_R = R_LSQ(1,NE)
     R12_R = R_LSQ(2,NE)
     R22_R = R_LSQ(3,NE)
     
     !Part 24:
     WX_R = DELX_R/(R11_R**2)-(R12_R/(R11_R*(R22_R**2)))*(DELY_R-DELX_R*(R12_R/R11_R))
     WY_R = (1.0/(R22_R**2))*(DELY_R-DELX_R*(R12_R/R11_R))
     

     !Part 25:
     GWNP1(1,1,NE)=GWNP1(1,1,NE) + WX_R*(R_L-R_R) 
     GWNP1(2,1,NE)=GWNP1(2,1,NE) + WY_R*(R_L-R_R)            
             
     GWNP1(1,2,NE)=GWNP1(1,2,NE) + WX_R*(U_L-U_R)
     GWNP1(2,2,NE)=GWNP1(2,2,NE) + WY_R*(U_L-U_R)             
               
     GWNP1(1,3,NE)=GWNP1(1,3,NE) + WX_R*(V_L-V_R)
     GWNP1(2,3,NE)=GWNP1(2,3,NE) + WY_R*(V_L-V_R)             
              
     GWNP1(1,4,NE)=GWNP1(1,4,NE) + WX_R*(P_L-P_R)
     GWNP1(2,4,NE)=GWNP1(2,4,NE) + WY_R*(P_L-P_R)            
     
 END DO  
 
 
!********************************************************************************************* 
 End
!###########################################################################################
