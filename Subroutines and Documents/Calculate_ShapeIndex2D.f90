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
!// Date: May., 15, 2016                                                                   //!
!// Developed by: M. Mohammadi, Mechanical Eng., Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Calculate_ShapeIndex2D(Dim,NC,X,Y,Corn,F_shape)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,X,Y,Corn
 Intent(Out  )::F_shape

 Integer::Dim,I,NC,P1,P2,P3
 Real(8)::X21,X31,Y21,Y31,Alpha,T
 Integer,Dimension(1:4,1:Dim)::Corn
 Real(8),Dimension(1:Dim)::X,Y,F_shape
!*********************************************************************************************
!Part 1: 
 DO I=1,NC
    P1 = Corn(1,I)
    P2 = Corn(2,I)
    P3 = Corn(3,I)

    X21=X(P2)-X(P1)
    X31=X(P3)-X(P1)
    Y21=Y(P2)-Y(P1)
    Y31=Y(P3)-Y(P1)
  
   !Part 3: 
    Alpha=SQRT((X21**2+Y21**2)*(X31**2+Y31**2)-(X21*X31+Y21*Y31)**2)

   !Part 4: 
    F_shape(I)=SQRT(3.)*Alpha/((X21**2+Y21**2)+(X31**2+Y31**2)-(X21*X31+Y21*Y31))

 ENDDO
!*********************************************************************************************
 End
!###########################################################################################