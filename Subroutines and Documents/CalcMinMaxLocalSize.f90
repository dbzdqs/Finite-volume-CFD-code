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
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine CalcMinMaxLocalSize(Metric,MinVal,MaxVal,EDir)
 Implicit None
!*********************************************************************************************
 Integer::I,EDir
 Real(8),Dimension(1:2,1:2)::Metric
 Real(8)::MinVal,MaxVal,TempReal1,TempReal2
 Real(8),Dimension(1:2,1:12)::Directions
 Logical::flag1
!*********************************************************************************************
 
!just 180 degree(1 to 6)
!Part 1:
 Directions(1,1) = 1.0
 Directions(2,1) = 0.0
 
 Directions(1,2) = 1.5/1.8027756377319946465596106337352
 Directions(2,2) = 1.0/1.8027756377319946465596106337352
 
 Directions(1,3) = 1.0/sqrt(5.0)
 Directions(2,3) = 2.0/sqrt(5.0)
 
 Directions(1,4) = 0.0
 Directions(2,4) = 1.0
 
 Directions(1,5) = -1.0/sqrt(5.0)
 Directions(2,5) = 2.0/sqrt(5.0)
 
 Directions(1,6) = -1.5/1.8027756377319946465596106337352
 Directions(2,6) = 1.0/1.8027756377319946465596106337352
 
 Directions(1,7) = -1.0
 Directions(2,7) = 0.0
 
 Directions(1,8) = -1.5/1.8027756377319946465596106337352
 Directions(2,8) = -1.0/1.8027756377319946465596106337352
 
 Directions(1,9) = -1.0/sqrt(5.0)
 Directions(2,9) = -2.0/sqrt(5.0)
 
 Directions(1,10) = 0.0
 Directions(2,10) = -1.0
 
 Directions(1,11) = 1.0/sqrt(5.0)
 Directions(2,11) = -2.0/sqrt(5.0)
 
 Directions(1,12) = 1.5/1.8027756377319946465596106337352
 Directions(2,12) = -1.0/1.8027756377319946465596106337352
 
!Part 2:
 flag1 = .TRUE.
 Do I=1,12
     
   !Part 3:
    TempReal1 = sqrt( (Directions(1,I)**2) + (Directions(2,I)**2) )   
    TempReal2 =  Directions(1,I)*( Directions(1,I)*Metric(1,1) + Directions(2,I)*Metric(1,2) )
    TempReal2 =  TempReal2 + Directions(2,I)*( Directions(1,I)*Metric(2,1) + Directions(2,I)*Metric(2,2) )
    TempReal1 = TempReal1/sqrt(TempReal2)
     
   !Part 4:
    If(MinVal>TempReal1 .OR. flag1) MinVal = TempReal1
    If(MaxVal<TempReal1 .OR. flag1)Then
     MaxVal = TempReal1
     EDir   = I
    Endif
     
    flag1 = .FALSE.
 EndDo
!*********************************************************************************************
 End
!###########################################################################################