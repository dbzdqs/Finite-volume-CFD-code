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
!// Date: Feb., 18, 2017                                                                   //!
!// Developed by: M. Mohammadi, Mechanical Eng., Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Calculate_SkewIndex2D(Dim,NC,X,Y,Corn,Skew)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,X,Y,Corn
 Intent(Out  )::Skew

 Integer::Dim,I,NC,P1,P2,P3
 Real(8):: Pi,U11,U12,U21,U22
 Integer,Dimension(1:4,1:Dim)::Corn
 Real(8),Dimension(1:Dim)::X,Y,Skew
 Real(8),Dimension(1:3):: Angle 
!*********************************************************************************************
!Part 1: 
 Pi=2*ASIN(1.)
 
 DO I=1,NC

    P1 = Corn(1,I)
    P2 = Corn(2,I)
    P3 = Corn(3,I)

   !Part 2: 
    U11=X(P2)-X(P1)
    U12=Y(P2)-Y(P1)
    U21=X(P3)-X(P1)
    U22=Y(P3)-Y(P1)

    Angle(1)=ACOS((U11*U21+U12*U22)/(SQRT(U11**2+U12**2)*SQRT(U21**2+U22**2)))

    U11=X(P1)-X(P2)
    U12=Y(P1)-Y(P2)
    U21=X(P3)-X(P2)
    U22=Y(P3)-Y(P2)

    Angle(2)=ACOS((U11*U21+U12*U22)/(SQRT(U11**2+U12**2)*SQRT(U21**2+U22**2)))
    
    U11=X(P1)-X(P3)
    U12=Y(P1)-Y(P3)
    U21=X(P2)-X(P3)
    U22=Y(P2)-Y(P3)

    Angle(3)=ACOS((U11*U21+U12*U22)/(SQRT(U11**2+U12**2)*SQRT(U21**2+U22**2)))

   !Part 3: 
    Skew(I)=max((MAXVAL(Angle)-Pi/3)/(Pi-Pi/3),(Pi/3-MINVAL(Angle))/(Pi/3))

 End Do
!*********************************************************************************************
 End
!###########################################################################################