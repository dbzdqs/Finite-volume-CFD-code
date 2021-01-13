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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ResidualSmoothing_expilicit(Dim,NC,NF1,NF2,IDS,Eps,Residual)
 Implicit None
!**********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,IDS,Eps
 Intent(InOut)::Residual

 Integer::Dim,I,ME,NE,NC,NF1,NF2
 Real(8)::Eps
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::Residual,sig
!**********************************************************************************************
!Part 1:
 sig(:,:)=0

!Part 2:
 Do I=NF1+1,NF2 

  !Part 3: 
   ME=IDS(1,I)
   NE=IDS(2,I)
   
  !Part 4: 
   sig(1,ME)=sig(1,ME)+Eps*(Residual(1,NE)-Residual(1,ME))
   sig(2,ME)=sig(2,ME)+Eps*(Residual(2,NE)-Residual(2,ME))
   sig(3,ME)=sig(3,ME)+Eps*(Residual(3,NE)-Residual(3,ME))
   sig(4,ME)=sig(4,ME)+Eps*(Residual(4,NE)-Residual(4,ME))
   sig(1,NE)=sig(1,NE)+Eps*(Residual(1,ME)-Residual(1,NE))
   sig(2,NE)=sig(2,NE)+Eps*(Residual(2,ME)-Residual(2,NE))
   sig(3,NE)=sig(3,NE)+Eps*(Residual(3,ME)-Residual(3,NE))
   sig(4,NE)=sig(4,NE)+Eps*(Residual(4,ME)-Residual(4,NE))
  
 End Do

!Part 5:
 Do I=1,NC
    Residual(1,I)=Residual(1,I)+sig(1,I)
    Residual(2,I)=Residual(2,I)+sig(2,I)
    Residual(3,I)=Residual(3,I)+sig(3,I)
    Residual(4,I)=Residual(4,I)+sig(4,I)
 End Do
!********************************************************************************************** 
 END 
!##############################################################################################
