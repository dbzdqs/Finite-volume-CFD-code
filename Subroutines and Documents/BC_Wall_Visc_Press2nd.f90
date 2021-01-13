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
 Subroutine BC_Wall_Visc_Press2nd(Dim,NFW1,NFW2,IDS,GM,P,WB,WNP1,X,Y,XC,YC,GWNP1,Limit)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,GM,P,WNP1,X,Y,XC,YC,GWNP1,Limit
 Intent(Out  )::WB

 Integer::Dim,I,NFW1,NFW2,ME
 Real(8)::GM,GM1,PB
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:5,1:Dim)::WB
 
 Real(8),Dimension(1:4,1:Dim)::WNP1
 
 Integer::P1,P2
 Real(8)::X_L,Y_L,X_P1,Y_P1,X_P2,Y_P2,XM_EDG,YM_EDG,DELX_L,DELY_L,RB
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:Dim)::X,Y,XC,YC
 
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8)::Phi
!*********************************************************************************************	
 !Part 1: 
 GM1= GM-1.

 !Part 2:
 DO I=NFW1+1,NFW2

    !Part 3:
    ME = IDS(1,I)
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
    Phi=Limit(1,ME)
    RB = WNP1(1,ME) + Phi * (GWNP1(1,1,ME)*DELX_L+GWNP1(2,1,ME)*DELY_L)
    
    Phi=Limit(4,ME)
    PB = P(ME) + Phi * (GWNP1(1,4,ME)*DELX_L+GWNP1(2,4,ME)*DELY_L)
    
	!Part 8:
    WB(1,I) = RB
    WB(2,I) = 0.0
    WB(3,I) = 0.0
    WB(4,I) = PB/GM1
    WB(5,I) = PB

 End do
!*********************************************************************************************
 End
!###########################################################################################
