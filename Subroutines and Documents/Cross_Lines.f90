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
 Subroutine Cross_Lines(X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros,Intersect)
 Implicit None
!*********************************************************************************************
 Intent(In   )::X1,X2,X3,X4,Y1,Y2,Y3,Y4
 Intent(Out  )::Xcros,Ycros,Intersect

 Real(8)::X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros,M12,M34,Eps,D1,D2,D12,D3,D4,D34
 Integer::Intersect
!*********************************************************************************************
!Part 1:
 Eps=0.00000001
 Intersect=-1
 If( Dabs(X2-X1)<Eps .and. Dabs(X4-X3)<Eps ) Goto 100

!Part 2:
 If( Dabs(X2-X1)<Eps )Then
  Xcros = X1
  M34=(Y4-Y3)/(X4-X3)
  Ycros = M34*Xcros -M34*X3 + Y3

!Part 3:
 Elseif( dabs(X4-X3)<Eps )Then
  Xcros = X3
  M12=(Y2-Y1)/(X2-X1)
  Ycros = M12*Xcros -M12*X1 + Y1

!Part 4:
 Else
  M12=(Y2-Y1)/(X2-X1)
  M34=(Y4-Y3)/(X4-X3)
 
  If( dabs(M34-M12)<Eps )M12=M34+Eps
  Xcros = (Y1-Y3+M34*X3-M12*X1)/(M34-M12)
  Ycros = M12*Xcros -M12*X1 + Y1
 Endif

!Part 5:
 D1 =Sqrt( (X1-Xcros)**2 + (Y1-Ycros)**2 )
 D2 =Sqrt( (X2-Xcros)**2 + (Y2-Ycros)**2 )
 D12=Sqrt( (X1-X2)**2 + (Y1-Y2)**2 )

!Part 6: 
 D3 =Sqrt( (X3-Xcros)**2 + (Y3-Ycros)**2 )
 D4 =Sqrt( (X4-Xcros)**2 + (Y4-Ycros)**2 )
 D34=Sqrt( (X3-X4)**2 + (Y3-Y4)**2 )

!Part 7: 
 If( Dabs(D34-D3-D4)<Eps .and. Dabs(D12-D1-D2)<Eps ) Intersect=1 
!*********************************************************************************************
100 End 
!###########################################################################################


   
