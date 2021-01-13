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
!// Date: Oct., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Kw_BC3D(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,MR,NX,NY,NZ,&
                    DW,Mu,WB,WNP1,WTNP1,WTB)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NX,NY,NZ,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,MR,DW,&
                 Mu,WB,WNP1,WTNP1
 Intent (Out  )::WTB

 Integer::Dim,I,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,ME
 Real(8)::MR,U,V,W,Q
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::Mu,DW,NX,NY,NZ
 Real(8),Dimension(1:2,1:Dim)::WTB,WTNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5,1:Dim)::WNP1
!*********************************************************************************************
!Part 1:
 DO I=NFI1+1,NFI2
 	WTB(1,I) = 1.0e-6
    WTB(2,I) = 5.0 
 END do
  
!Part 2:  
 DO I=NFO1+1,NFO2
    ME       = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 3:
 DO I=NFW1+1,NFW2
    ME       = IDS(1,I)
 	WTB(1,I) = 0.0
    WTB(2,I) = MR * WB(1,I)*(60.0*Mu(ME)/ (WNP1(1,ME)*0.075*DW(ME)*DW(ME)) ) 
 END do
  
!Part 4:
 DO I=NFS1+1,NFS2
    ME       = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 5:
 DO I=NFF1+1,NFF2
    ME  = IDS(1,I)
    
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)
    
    Q  = U*NX(I)+V*NY(I)+W*NZ(I)
        
    IF( Q<=0. )Then  
 	 WTB(1,I) = 1.0e-6
     WTB(2,I) = 5.0
    Else
 	 WTB(1,I) = WTNP1(1,ME)
     WTB(2,I) = WTNP1(2,ME)    
    End IF
 END Do
!*********************************************************************************************
 End
!###########################################################################################

