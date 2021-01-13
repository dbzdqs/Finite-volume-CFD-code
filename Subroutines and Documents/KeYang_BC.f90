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
 Subroutine KeYang_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,DW,Wb,Wnp1,WTNP1,IDS,MR,Mu,WTB)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,DW,Wb,Wnp1,WTNP1,IDS,MR,Mu
 Intent(Out  )::WTB

 Integer::Dim,I,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,ME,P1,P2
 Real(8)::MR,U,V,Q
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::Mu,DW,NX,NY
 Real(8),Dimension(1:2,1:Dim)::WTB,WTNP1
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:4,1:Dim)::Wnp1
!**********************************************************************************************	
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
    ME  = IDS(1,I)
 	WTB(1,I) = 0.0
    WTB(2,I) = 2.0*MR*Mu(ME)*(dsqrt(dabs(WTNP1(1,ME)))/DW(ME))*(dsqrt(dabs(WTNP1(1,ME)))/DW(ME))
  END do
  
  !Part 4:
  DO I=NFS1+1,NFS2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
  END do
  
  !part 5:
  DO I=NFF1+1,NFF2
    ME  = IDS(1,I)
    
    U = Wb(2,I)/Wb(1,I)
    V = Wb(3,I)/Wb(1,I)
    
    Q  = U*NX(I)+V*NY(I)
        
    if (Q<=0.0) then
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