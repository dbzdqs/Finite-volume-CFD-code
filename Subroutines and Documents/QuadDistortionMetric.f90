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
Function QuadDistortionMetric(Dim,Corn,X,Y,Elm,COINCIDENT_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Elm,COINCIDENT_TOLERANCE

Integer::Dim,Elm,A,B,C,D,I,J,negval,InvertedTris
Integer,Dimension(1:Dim,1:4)::Corn
Logical::hasAngleLessThanSixDegrees,hasTwoCoincidentNode
Real(8)::TriDistortionMetric,QuadDistortionMetric,Alpha1,Alpha2,Alpha3,Alpha4,COINCIDENT_TOLERANCE
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., J. R. Tristano and M. L. Staten (1998), "An approach to combined Laplacian and
!--->>> optimization-based smoothing for triangular, quadrilateral and tetrahedral meshes"

InvertedTris = 0

A = Corn(Elm,1)
B = Corn(Elm,2)
C = Corn(Elm,3)
D = Corn(Elm,4)

!Part 1:

!-------------------->>> Notice: Suppose ABCD is the Quad then: <<<-------------------------

Alpha1 = TriDistortionMetric(Dim,X,Y,A,B,C)
Alpha2 = TriDistortionMetric(Dim,X,Y,A,C,D)
Alpha3 = TriDistortionMetric(Dim,X,Y,B,C,D)
Alpha4 = TriDistortionMetric(Dim,X,Y,A,B,D)

if(Alpha1 < 0) InvertedTris = InvertedTris + 1

if(Alpha2 < 0) InvertedTris = InvertedTris + 1

if(Alpha3 < 0) InvertedTris = InvertedTris + 1

if(Alpha4 < 0) InvertedTris = InvertedTris + 1
 
!Part 2:

if(hasAngleLessThanSixDegrees(Dim,Corn,X,Y,Elm) .Or. hasTwoCoincidentNode(Dim,Corn,X,Y,Elm,COINCIDENT_TOLERANCE) .Or. InvertedTris == 2) then

	negval = 1

elseif(InvertedTris == 3) then

	negval = 2

elseif(InvertedTris == 4) then

	negval = 3

else

	negval = 0

endif

QuadDistortionMetric = DMIN1(Alpha1,Alpha2,Alpha3,Alpha4) - negval

!===========================================================================================
End Function QuadDistortionMetric
!*********************************************************************************************
