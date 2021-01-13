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
!// Supervisor:                                                                            //!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Date: Feb., 10, 2018                                                                   //!
!// Supervisor:                                                                            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GenGosianPoint(Dim,NF,IDS,X,Y,XF,YF)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,IDS,X,Y
 Intent(Out  )::XF,YF

 Integer::Dim,I,P1,P2,NF
 Real(8)::Coeff,Temp,TempX,TempY
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:2,1:Dim)::XF,YF
!*********************************************************************************************	
!Part1:
 Coeff = Sqrt(1.0 / 3.0)
 
!Part 2:
 DO I=1,NF

  !Part 3:
   P1 = IDS(3,I)
   P2 = IDS(4,I)
   
  !Part 4:
   TempX = ( X(P2) + X(P1) ) / 2.0
   TempY = ( Y(P2) + Y(P1) ) / 2.0
   
  !Part 5:    
   XF(1,I) = (1.0-Coeff) * TempX  + Coeff *  X(P1)
   YF(1,I) = (1.0-Coeff) * TempY  + Coeff *  Y(P1)
   
  !Part 6:
   XF(2,I) = (1.0-Coeff) * TempX  + Coeff *  X(P2)
   YF(2,I) = (1.0-Coeff) * TempY  + Coeff *  Y(P2)

 End Do
!*********************************************************************************************
 End
!###########################################################################################