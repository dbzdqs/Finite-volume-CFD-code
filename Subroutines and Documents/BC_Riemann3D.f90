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
 Subroutine BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,Ro0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
 Implicit None
!**********************************************************************************************
 Intent(In   )::Dim,NFF1,NFF2,GM,U0,V0,W0,P0,Ro0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P
 Intent(Out  )::WB

 Integer::Dim,I,NFF1,NFF2,ME
 Real(8)::GM,GM1,DS,Nxx,Nyy,Nzz,UI,VI,WI,CI,MI,QNI,RI,R0,SI,QN0,S0,QN,S,UB,VB,WBB,SB,ROB,PB1,U0,V0,W0,P0,Ro0,C0
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::Nx,Ny,Nz,DA,P
 Real(8),Dimension(1:6,1:Dim)::WB
!**********************************************************************************************	
 GM1= GM-1.
 
!Part 1:
 DO I=NFF1+1,NFF2

   !Part 2:
    DS  = DA(I)
	Nxx = Nx(I)/DS    
	Nyy = Ny(I)/DS    
	Nzz = Nz(I)/DS

   !Part 3:
	ME  = IDS(1,I)
	UI  = WNP1(2,ME)/WNP1(1,ME)
	VI  = WNP1(3,ME)/WNP1(1,ME)
	WI  = WNP1(4,ME)/WNP1(1,ME)

	CI  = DSQRT(DABS(GM*P(ME)/WNP1(1,ME)))
	MI  = DSQRT(UI*UI+VI*VI+WI*WI)/CI

   !Part 4:
	QNI = UI*Nxx+VI*Nyy+WI*Nzz
	RI  = QNI+2*CI/GM1
	SI  = P(ME)/WNP1(1,ME)**GM

   !Part 5:
	QN0 = U0*Nxx+V0*Nyy+W0*Nzz
	R0  = QN0-2*C0/GM1
	S0  = P0/Ro0**GM

   !Part 6:
	QN = 0.5 *(RI+R0)
	S  = 0.25*GM1*(RI-R0)
	
   !Part 7:
	IF(MI >= 1.0)THEN       !SUPERSONIC BOUNDARY CONDITION

    !Part 8:
	 IF( QN<0.0 )THEN        !INPUT BOUNDARY
	  UB  = U0
	  VB  = V0
	  WBB = W0
	  ROB = Ro0 
	  PB1 = P0
      
    !Part 9:
	 ELSE                  !OUTPUT BOUNDARY
	  UB  = UI
	  VB  = VI
	  WBB = WI  
	  ROB = WNP1(1,ME) 
	  PB1 = P(ME)
	 END IF

   !Part 10:
    ELSEIF(MI < 1.0)THEN  !SUBSONIC BOUNDARY CONDITION
         
    !Part 11:
	 IF( QN<0.0 )THEN          !INPUT BOUNDARY
	  UB  = U0+(QN-QN0)*Nxx
	  VB  = V0+(QN-QN0)*Nyy
	  WBB = W0+(QN-QN0)*Nzz
	  SB  = S0
	  ROB = (S*S/GM/SB)**(1/GM1)
	  PB1 = S*S*ROB/GM
      
    !Part 12:
	 ELSE                      !OUTPUT BOUNDARY
	  UB  = UI+(QN-QNI)*Nxx
	  VB  = VI+(QN-QNI)*Nyy
	  WBB = WI+(QN-QNI)*Nzz
	  SB  = SI
	  ROB = (S*S/GM/SB)**(1/GM1)
	  PB1 = S*S*ROB/GM
	 END IF

	END IF

   !Part 13:
	WB(1,I) =	ROB
	WB(2,I) =	ROB*UB
	WB(3,I) =	ROB*VB
	WB(4,I) =	ROB*WBB
	WB(5,I) =	PB1/(GM1) + 0.5*ROB*(UB*UB + VB*VB + WBB*WBB)
	WB(6,I) =	PB1

 END DO
!**********************************************************************************************
 End
!##############################################################################################


