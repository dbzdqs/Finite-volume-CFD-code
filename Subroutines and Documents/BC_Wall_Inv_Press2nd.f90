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
!// Date: Feb., 15, 2017                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine	 BC_Wall_Inv_Press2nd(Dim,NFW1,NFW2,NX,NY,DA,IDS,GM,WNP1,P,WB,X,Y,XC,YC,GWNP1,Limit)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NFW1,NFW2,GM,IDS,Wnp1,NX,NY,DA,P,X,Y,XC,YC,GWNP1,Limit
 Intent (Out  )::Wb

 Integer::Dim,I,NFW1,NFW2,ME,P1,P2
 Real(8)::GM,GM1,U,V,CC,MLoc,C0,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,STH,CTH,QN0,QT0,RI0,PE,RE,CE,&
          VE,QNE,QTE,QNB,QTB,RIE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0,ITT,ITP,AII,&
		  HTE,RIB,HTB,EP1,EP2,EP3,C1,C2,MB,TB,TE
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:5,1:Dim)::Wb
 
 
 Real(8)::X_L,Y_L,X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG,DELX_L,DELY_L
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y,XC,YC
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8)::Phi
!*********************************************************************************************
 !Part 1:
 GM1= GM-1
 
 !Part 2:
 DO I=NFW1+1,NFW2

    !Part 3:
    ME  = IDS(1,I)
    X_L=XC(ME)
    Y_L=YC(ME)
    
    !Part 4:
    P1=IDS(3,I)          
    P2=IDS(4,I)          
                        
    X_P1=X(P1)
    Y_P1=Y(P1)
    
    X_P2=X(P2)
    Y_P2=Y(P2)
    
    !Part 5:
    XM_EDG=0.5*(X_P1+X_P2)
    YM_EDG=0.5*(Y_P1+Y_P2)
    
    !Part 6:
    DELX_L = (XM_EDG-X_L)
    DELY_L = (YM_EDG-Y_L)
    
    !Part 7:
    STH = NX(I)/DA(I)
    CTH = NY(I)/DA(I)

    !Part 8:
    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
    
    !Part 9:
	QNE = UE*STH+VE*CTH
    QTE =-UE*CTH+VE*STH
    
    !Part 10:
    QNB=0
	QTB=QTE
    
    !Part 11:
	UB=STH*QNB-CTH*QTB
	VB=CTH*QNB+STH*QTB
	
	!Part 12:
    Phi=Limit(1,ME)
    RB = WNP1(1,ME) + Phi * (GWNP1(1,1,ME)*DELX_L+GWNP1(2,1,ME)*DELY_L)
    
    Phi=Limit(4,ME)
    PB = P(ME) + Phi * (GWNP1(1,4,ME)*DELX_L+GWNP1(2,4,ME)*DELY_L)
    
    !Part 13:
	Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = PB/GM1 + 0.5*RB*(UB*UB + VB*VB)
    Wb(5,I) = PB

	END do

!*********************************************************************************************
 End
!###########################################################################################