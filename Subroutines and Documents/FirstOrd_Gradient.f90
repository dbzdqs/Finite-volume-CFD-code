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
 Subroutine FirstOrd_Gradient(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P
 Intent(Out  )::GWNP1

 Integer::Dim,I,J,NC,NF1,NF2,NF,ME,NE
 Real(8)::NXX,NYY,R,U,V,Pr,R_NX,U_NX,V_NX,P_NX,R_NY,U_NY,V_NY,P_NY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,A  ,P    
 Real(8),Dimension(1:4,1:Dim)::WNP1      
 Real(8),Dimension(1:5,1:Dim)::WB    
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
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
    R_NX = R*NXX
    U_NX = U*NXX
    V_NX = V*NXX
    P_NX = Pr*NXX

    R_NY = R*NYY
    U_NY = U*NYY
    V_NY = V*NYY
    P_NY = Pr*NYY

    GWNP1(1,1,ME)=GWNP1(1,1,ME) + R_NX
    GWNP1(2,1,ME)=GWNP1(2,1,ME) + R_NY

    GWNP1(1,2,ME)=GWNP1(1,2,ME) + U_NX
    GWNP1(2,2,ME)=GWNP1(2,2,ME) + U_NY

    GWNP1(1,3,ME)=GWNP1(1,3,ME) + V_NX
    GWNP1(2,3,ME)=GWNP1(2,3,ME) + V_NY
      
    GWNP1(1,4,ME)=GWNP1(1,4,ME) + P_NX
    GWNP1(2,4,ME)=GWNP1(2,4,ME) + P_NY

 End Do

!Part 6: 
 DO I=NF1+1,NF2

   !Part 7: 
    ME = IDS(1,I)
    NE = IDS(2,I)
    
	NXX = NX(I)
    NYY = NY(I)

   !Part 8:
	R = 0.5*( WNP1(1,ME)            + WNP1(1,NE)            )
    U = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )
    V = 0.5*( WNP1(3,ME)/WNP1(1,ME) + WNP1(3,NE)/WNP1(1,NE) )
    Pr= 0.5*( P(ME)                 + P(NE)                 )

    R_NX = R*NXX
    U_NX = U*NXX
    V_NX = V*NXX
    P_NX = Pr*NXX

    R_NY = R*NYY
    U_NY = U*NYY
    V_NY = V*NYY
    P_NY = Pr*NYY

    GWNP1(1,1,ME)=GWNP1(1,1,ME) + R_NX
    GWNP1(2,1,ME)=GWNP1(2,1,ME) + R_NY             
    GWNP1(1,1,NE)=GWNP1(1,1,NE) - R_NX
    GWNP1(2,1,NE)=GWNP1(2,1,NE) - R_NY
   
    GWNP1(1,2,ME)=GWNP1(1,2,ME) + U_NX
    GWNP1(2,2,ME)=GWNP1(2,2,ME) + U_NY             
    GWNP1(1,2,NE)=GWNP1(1,2,NE) - U_NX
    GWNP1(2,2,NE)=GWNP1(2,2,NE) - U_NY
	  
    GWNP1(1,3,ME)=GWNP1(1,3,ME) + V_NX
    GWNP1(2,3,ME)=GWNP1(2,3,ME) + V_NY             
    GWNP1(1,3,NE)=GWNP1(1,3,NE) - V_NX
    GWNP1(2,3,NE)=GWNP1(2,3,NE) - V_NY	   

    GWNP1(1,4,ME)=GWNP1(1,4,ME) + P_NX
    GWNP1(2,4,ME)=GWNP1(2,4,ME) + P_NY             
    GWNP1(1,4,NE)=GWNP1(1,4,NE) - P_NX
    GWNP1(2,4,NE)=GWNP1(2,4,NE) - P_NY
 
 End Do

!Part 9: 
 DO I=1,NC
    Do J=1,4    
       GWNP1(1,J,I)=GWNP1(1,J,I)/A(I)    
       GWNP1(2,J,I)=GWNP1(2,J,I)/A(I)      
    End Do
 End Do 
!********************************************************************************************* 
 End
!###########################################################################################
