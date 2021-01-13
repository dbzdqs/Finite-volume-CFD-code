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
SUBROUTINE Write_MeshQuality2D(Dim,NC,Skew,Aspect,Smoothness,F_size,F_shape,F_compound)
 IMPLICIT NONE
!*********************************************************************************************
 INTEGER::Dim,NC,I
 Real*8,DIMENSION(1:Dim)::Skew,Aspect,F_size,F_shape,F_compound,Smoothness
!*********************************************************************************************
!Part 1: 
 OPEN(1,File='Skew.Plt')
 OPEN(2,File='Aspect.Plt')
 OPEN(3,File='Smoothness.Plt')
 OPEN(4,File='size_index.Plt') 
 OPEN(5,File='shape_index.Plt')
 OPEN(6,File='compound_index.Plt')

!Part 2: 
 DO I=1,NC
	Write(1,*) Skew(I)
 END DO

 DO I=1,NC
	Write(2,*) Aspect(I)
 END DO

 DO I=1,NC
	Write(3,*) Smoothness(I)
 END DO

 DO I=1,NC
	Write(4,*) F_size(I)
 END DO

 DO I=1,NC
	Write(5,*) F_shape(I)
 END DO

 DO I=1,NC
	Write(6,*) F_compound(I)
 END DO
 
 CLOSE(1)
 CLOSE(2)
 CLOSE(3)
 CLOSE(4)
 CLOSE(5)
 CLOSE(6)
!*********************************************************************************************
 END