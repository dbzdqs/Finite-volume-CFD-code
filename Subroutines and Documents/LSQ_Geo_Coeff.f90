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
 Subroutine LSQ_GEO_COEFF(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,XC,YC,X,Y,R_LSQ)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,XC,YC,X,Y
 Intent(Out  )::R_LSQ

 Integer::Dim,I,J,NC,NF1,NF2,NF,ME,NE
 Real(8)::NXX,NYY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,A        
 
 Real(8),Dimension(1:Dim)::X,Y,XC,YC
 Real(8),Dimension(1:3,1:Dim)::R_LSQ
 Real(8)::X_L,Y_L,X_R,Y_R,X_P1,X_P2,Y_P1,Y_P2,XM_EDG,YM_EDG,DELX_L,DELY_L,DELX_R,DELY_R
 
 Real(8)::R11_L,R12_L,R22_L,R11_R,R12_R,R22_R,WX_L,WX_R,WY_L,WY_R
 Integer::P1,P2

 Real(8)::XP_L,YP_L,SLOPE,CNST

!*********************************************************************************************	
!Part 1:
 R_LSQ=0.0
 	
!Part 2: 
 DO I=NF2+1,NF   

   !Part 3: 
    ME = IDS(1,I) 
   
    !Part 4:
    X_L=XC(ME)
    Y_L=YC(ME)
    
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                         
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    !Part 5:     
    IF ((X_P2-X_P1)/=0.0 .AND. (Y_P2-Y_P1)/=0) THEN
        
        SLOPE = (Y_P2-Y_P1)/(X_P2-X_P1)
        CNST = Y_P1-SLOPE*X_P1
        
        XP_L = (2.0*(Y_L-CNST)+X_L*((1.0/SLOPE)-SLOPE))/(SLOPE+(1.0/SLOPE))
        YP_L = Y_L-(1.0/SLOPE)*(XP_L-X_L)
        
    !Part 6:    
    ELSE IF ((X_P2-X_P1)==0.0) THEN
        
        XP_L = 2.0*X_P1-X_L
        YP_L = Y_L
   
    !Part 7:
    ELSE IF ((Y_P2-Y_P1)==0.0) THEN
        
        XP_L = X_L
        YP_L = 2.0*Y_P1-Y_L
        
    END IF
    
    !Part 8:
    DELX_L = (XP_L-X_L)
    DELY_L = (YP_L-Y_L)
    
    !Part 9:
    R_LSQ(1,ME) = R_LSQ(1,ME)+DELX_L**2
    R_LSQ(2,ME) = R_LSQ(2,ME)+DELX_L*DELY_L
    R_LSQ(3,ME) = R_LSQ(3,ME)+DELY_L**2
    
    
 End Do

!Part 10: 
 DO I=NF1+1,NF2

   !Part 11: 
    ME = IDS(1,I)
    NE = IDS(2,I)
    
    !Part 12:
	X_L=XC(ME)
    Y_L=YC(ME)
    
    X_R=XC(NE)
    Y_R=YC(NE)
    
    !Part 13:
    DELX_L = X_R-X_L
    DELY_L = Y_R-Y_L
        
    DELX_R = X_L-X_R
    DELY_R = Y_L-Y_R
    
    !Part 14:
    R_LSQ(1,ME) = R_LSQ(1,ME)+DELX_L**2
    R_LSQ(2,ME) = R_LSQ(2,ME)+DELX_L*DELY_L
    R_LSQ(3,ME) = R_LSQ(3,ME)+DELY_L**2
    
    R_LSQ(1,NE) = R_LSQ(1,NE)+DELX_R**2
    R_LSQ(2,NE) = R_LSQ(2,NE)+DELX_R*DELY_R
    R_LSQ(3,NE) = R_LSQ(3,NE)+DELY_R**2
     
 END DO
 
    
!Part 15:    
DO I=1,NC
     
     R_LSQ(1,I) = DSQRT(R_LSQ(1,I))
     R_LSQ(2,I) = R_LSQ(2,I)/R_LSQ(1,I)
     R_LSQ(3,I) = DSQRT(R_LSQ(3,I)-R_LSQ(2,I)**2)
     
 END DO   
!********************************************************************************************* 
 End
!###########################################################################################
