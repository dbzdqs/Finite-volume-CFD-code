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
 Subroutine LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Phi,Vol,Phihat)
 Implicit None
!********************************************************************************************* 
 Intent(In ):: Dim,NC,NF1,NF2,IDS,Phi,Vol
 Intent(Out )::Phihat

 Integer::Dim,I,NC,NF1,NF2,ME,NE
 Real(8),Dimension(1:Dim)::Phi,Vol,Phihat
 Real(8),Dimension(1:Dim)::Num,Denom
 Integer,Dimension(1:6,1:Dim)::IDS
!********************************************************************************************* 
!Part 1:
 Do I=1,NC
    Num(I)=Phi(I)*Vol(I)
    Denom(I)=Vol(I)    
 End Do

!Part 2:
 Do I=NF1+1,NF2
    
   !Part 3:
    ME=IDS(1,I)
    NE=IDS(2,I)
    
   !Part 4:
    Num(ME)   = Num(ME)   + Phi(NE)*Vol(NE)
    Denom(ME) = Denom(ME) + Vol(NE) 
	   
    Num(NE)   = Num(NE)   + Phi(ME)*Vol(ME)
    Denom(NE) = Denom(NE) + Vol(ME)

 End Do

!Part 5:
 Do I=1,NC
    Phihat(I)=Num(I)/Denom(I)
 End Do
!*********************************************************************************************
 End
!###########################################################################################