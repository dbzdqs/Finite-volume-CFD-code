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
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeChien_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS
 Intent(Out  )::WTB

 Real(8)::U,V,Q
 Integer::Dim,I,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,ME,P1,P2
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY
 Real(8),Dimension(1:2,1:Dim)::WTB,WTNP1
 Real(8),Dimension(1:5,1:Dim)::WB
!*********************************************************************************************
!Part 1:
 DO I=NFI1+1,NFI2
 	WTB(1,I) = 1.0e-6
    WTB(2,I) = 4.5e-7 
 END do
  
!Part 2:  
 DO I=NFO1+1,NFO2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 3:
 DO I=NFW1+1,NFW2
 	WTB(1,I) = 0.0
    WTB(2,I) = 0.0 
 END do
  
!Part 4:
 DO I=NFS1+1,NFS2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 5:
 DO I=NFF1+1,NFF2
    ME  = IDS(1,I)
    
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    
    Q  = U*NX(I)+V*NY(I)
        
    if (Q>=0.0) then
 	WTB(1,I) = 1.0e-6
    WTB(2,I) = 4.5e-7
    else
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)    
    end if

 END do
!*********************************************************************************************
 End
!###########################################################################################
