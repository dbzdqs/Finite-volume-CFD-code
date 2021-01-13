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
!// Date: Nov., 15, 2014                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function TriDistortionMetric(Dim,X,Y,A,B,C)
Implicit None
!===========================================================================================
Intent(In)::A,B,C

Integer::Dim,A,B,C,I
Real(8)::temp1,temp2,TriDistortionMetric,ABx,ABy,CAx,CAy,CBx,CBy
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., J. R. Tristano and M. L. Staten (1998), "An approach to combined Laplacian and
!--->>> optimization-based smoothing for triangular, quadrilateral and tetrahedral meshes"
!Part 1:
!---------------------->>> Notice: Suppose ABC is the Triangle then: <<<--------------------

ABx = X(B) - X(A) !--------------------- x_coordinate of AB vector -------------------------
ABy = Y(B) - Y(A) !--------------------- y_coordinate of AB vector -------------------------

CAx = X(A) - X(C) !--------------------- x_coordinate of CA vector -------------------------
CAy = Y(A) - Y(C) !--------------------- y_coordinate of CA vector -------------------------

CBx = X(B) - X(C) !--------------------- x_coordinate of CB vector -------------------------
CBy = Y(B) - Y(C) !--------------------- y_coordinate of CB vector -------------------------

!------------------------ Calculating Surface Normal of Triangle --------------------------- 
 
temp1 = CAx*CBy - CBx*CAy !-- temp1 is the norm of 'Cross Product' of CA x CB
temp2 = (CAx*CAx + CAy*CAy) + (ABx*ABx + ABy*ABy) + (CBx*CBx + CBy*CBy)

if(temp1 < 0) then

	I = -1

elseif(temp1 > 0) then

	I = 1

endif

TriDistortionMetric = 2*SQRT(3.0)*(DABS(temp1)/temp2)*I

!===========================================================================================
End Function TriDistortionMetric
!*********************************************************************************************
