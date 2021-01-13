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
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function TriangleQuality(Dim,X,Y,A,B,C) !----->> Implementation Based on Lee (1994) <<------
Implicit None
!===========================================================================================
Intent(In)::A,B,C

Integer::Dim,A,B,C
Real(8)::temp1,temp2,TriangleQuality
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
!---------------------->>> Notice: Suppose ABC is the Triangle then: <<<--------------------
 
temp1 = ((X(A) - X(C))*(Y(B) - Y(C))) - ((X(B) - X(C))*(Y(A) - Y(C))) !-- temp1 is the norm of 'Cross Product' of CA x CB
temp2 = ((X(A) - X(C))*(X(A) - X(C)) + (Y(A) - Y(C))*(Y(A) - Y(C))) + ((X(B) - X(A))*(X(B) - X(A)) + (Y(B) - Y(A))*(Y(B) - Y(A))) + ((X(C) - X(B))*(X(C) - X(B)) + (Y(C) - Y(B))*(Y(C) - Y(B)))

TriangleQuality = 2*SQRT(3.0)*(temp1/temp2) 
!===========================================================================================
End Function TriangleQuality
!*********************************************************************************************
