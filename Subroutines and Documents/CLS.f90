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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine CLS(Dim,Corn,Neib,X,Y,Elm,V,COINCIDENT_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,Elm,V,COINCIDENT_TOLERANCE
Intent(InOut)::X,Y

Integer::Dim,V,Elm,E,A,B,C,NPC,TEC,QEC,I,J
Integer,Dimension(1:1000)::QElms,TElms,NPList
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::MoveIsAcceptable
Real(8)::COINCIDENT_TOLERANCE,x_value,y_value,Px,Py,TriDistortionMetric,QuadDistortionMetric,x_old,y_old
Real(8),Dimension(1:Dim)::X,Y
Real(8),Dimension(1:100)::PreTMetrics,PostTMetrics,PreQMetrics,PostQMetrics
!===========================================================================================

!--->>> Reference: Canann, S. A., J. R. Tristano and M. L. Staten (1998), "An approach to combined Laplacian and
!--->>> optimization-based smoothing for triangular, quadrilateral and tetrahedral meshes"

x_old = X(V) !-------------------- Initial V's x-coordinates -------------------------------
y_old = Y(V) !-------------------- Initial V's y-coordinates -------------------------------
!Part 1:
Call GetSurroundingElements(Dim,Corn,Neib,V,Elm,TElms,QElms,TEC,QEC)
Call GetNeibouringPoints(Dim,Corn,Neib,QElms,QEC,NPList,NPC,V,Elm)

!--------------------------- Calculate Initial Distortion Metrics --------------------------
!Part 2:
do J=1,QEC

	PreQMetrics(J) = QuadDistortionMetric(Dim,Corn,X,Y,QElms(J),COINCIDENT_TOLERANCE) 

end do

do J=1,TEC

    E = TElms(J)
    
    PreTMetrics(J) = TriDistortionMetric(Dim,X,Y,Corn(E,1),Corn(E,2),Corn(E,3))
    
end do

!Part 3:

!----------- Computing the location where Laplacian smoothing would place the node ---------

Call CenterOfPolygon(Dim,NPC,NPList,X,Y,x_value,y_value)

!Part 4:

!--------- Defining a Vector(P) from initial V coordinates to its new coordinates ----------

Px = x_value - x_old
Py = y_value - y_old

!-------------------------------- Moving V to new Coordinate -------------------------------

X(V) = x_old + Px 
Y(V) = y_old + Py 

!Part 5:

do I=1,20    
    
	!------------------ Compute distortion metric for neighboring elements -----------------

	do J=1,QEC

		PostQMetrics(J) = QuadDistortionMetric(Dim,Corn,X,Y,QElms(J),COINCIDENT_TOLERANCE) 

    end do
    
    do J=1,TEC

        E = TElms(J)
    
        PostTMetrics(J) = TriDistortionMetric(Dim,X,Y,Corn(E,1),Corn(E,2),Corn(E,3))
    
    end do

	!------- If the new location is 'acceptable' then break out of the loop ----------------
       
	if(MoveIsAcceptable(Dim,Corn,Neib,X,Y,TElms,TEC,QElms,QEC,PreTMetrics,PostTMetrics,PreQMetrics,PostQMetrics)) then

		exit

	else !--- If is not acceptable, then cut the proposed move distance in half and set this as the new location

		Px = Px/2
		Py = Py/2

		X(V) = x_old + Px 
		Y(V) = y_old + Py 

	endif

end do

!===========================================================================================
End Subroutine CLS
!*********************************************************************************************
