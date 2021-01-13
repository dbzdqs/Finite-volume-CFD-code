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
!// Date: Aug., 30, 2015                                                                   //!
!// Developed by: M. Mohammadi, Mechanical Eng., Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Calculate_AspectRatioIndex2D(Dim,NC,X,Y,Corn,AspectRatio)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,X,Y,Corn
 Intent(Out  )::AspectRatio

 Integer::Dim,I,NC,P1,P2,P3
 Real(8):: X21,X31,X32,Y21,Y31,Y32
 Integer,Dimension(1:4,1:Dim)::Corn
 Real(8),Dimension(1:Dim)::X,Y,AspectRatio
 Real(8),DIMENSION(1:3)::Length
!*********************************************************************************************
!Part 1: 
 DO I=1,NC

    P1 = Corn(1,I)
    P2 = Corn(2,I)
    P3 = Corn(3,I)

    X21=X(P2)-X(P1)
    X31=X(P3)-X(P1)
    X32=X(P3)-X(P1)
    Y21=Y(P2)-Y(P1)
    Y31=Y(P3)-Y(P1)
    Y32=Y(P3)-Y(P1)

   !Part 2: 
    Length(1) = SQRT( Y21*Y21 + X21*X21 )
    Length(2) = SQRT( Y31*Y31 + X31*X31 )
    Length(3) = SQRT( Y32*Y32 + X32*X32 )
   
   !Part 3: 
    AspectRatio(I)=MAXVAL(Length)/MINVAL(Length)

 End Do
!*********************************************************************************************
 End
!###########################################################################################