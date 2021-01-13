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
 Subroutine Quality_AreaLen(Dim,P1,P2,P3,X,Y,Qual)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,P1,P2,P3,X,Y
 Intent(Out  )::Qual

 Real(8)::Edge_L,Qual,D12,D23,D31,Area,X12,Y12,X23,Y23,X13,Y13
 Integer::Dim,P1,P2,P3
 Real(8),Dimension(1:Dim)::X,Y
!*********************************************************************************************
!Part 1:
 X12 = X(P2)-X(P1) ; Y12 = Y(P2)-Y(P1)
 X13 = X(P3)-X(P1) ; Y13 = Y(P3)-Y(P1)  
 X23 = X(P3)-X(P2) ; Y23 = Y(P3)-Y(P2) 
 
!Part 2:
 Area =  X12*Y13 - Y12*X13

!Part 3:
 D12 =  X12*X12 + Y12*Y12 
 D23 =  X23*X23 + Y23*Y23 
 D31 =  X13*X13 + Y13*Y13 

!Part 4:
 Edge_L = ( D12 + D23 + D31 )**1.5

!Part 5:
 Qual = 6*Sqrt(2.)*( Area / Edge_L )
!*********************************************************************************************
 End
!###########################################################################################





