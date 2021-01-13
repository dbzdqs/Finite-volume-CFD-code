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
!// Date: Mar., 05, 2013                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BoundPointLabeling(Dim,IDS,NR,NFR,BC,NBP,IBP)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,IDS,NR,NFR,BC
 Intent(Out  )::NBP,IBP

 Integer::DIM,NBP,NR,NBE,I,J,P
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:Dim)::IBP
!*********************************************************************************************
!Part 1:
 NBP = 0
 NBE = 0

!Part 2:
 Do I=1,NR
     
   !Part 3:
    if( BC(I)==1 ) goto 10   
    Do J=NBE+1,NBE+NFR(I)
        
      !Part 4:
       P        = IDS(3,J)
       NBP      = NBP + 1
	   IBP(NBP) = P
       
    End do
10  NBE = NBE + NFR(I)
 END DO
!*********************************************************************************************
 End
!###########################################################################################
