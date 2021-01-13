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
Subroutine AddNewPoint(tri_up,tri_down,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NK,p1,p2,NM,Nn) ! Tested
Implicit None
!===========================================================================================
Intent(In )::tri_up,tri_down,x_value,y_value,Dim,NK,p1,p2,NM,Fronts,States
Intent(Out)::X,Y
Intent(Inout)::NC,NP,Neib,Corn,FrontEdges,Nn

Integer::Dim,Fronts,tri_up,tri_down,NK,NM,Nn,NC,NP,m1,m2,m3,m4,p1,p2,N1,N2,N3,N4,I,I1,I2,I3,J1,J2,J3
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Real(8)::x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
do I=1,3  ! Finding Neibour elements of two splitting triangles in the following loop

	if(Corn(tri_up,I) == P2) then
		N1 = Neib(tri_up,I)
	elseif(Corn(tri_up,I) == P1) then
		N2 = Neib(tri_up,I)
	endif

	if(Corn(tri_down,I) == P2) then
		N3 = Neib(tri_down,I)
	elseif(Corn(tri_down,I) == P1) then
		N4 = Neib(tri_down,I)
	endif

end do
!Part 2:
do I=1,3
	if(Corn(tri_up,I) == NM) I1 = I
	if(Corn(tri_up,I) == P1) I2 = I
	if(Corn(tri_up,I) == P2) I3 = I

	if(Corn(tri_down,I) == NK) J1 = I
	if(Corn(tri_down,I) == P2) J2 = I
	if(Corn(tri_down,I) == P1) J3 = I
end do

!Part 3:

m1 = tri_up
NC = NC+1
m2 = NC
m3 = tri_down
NC = NC +1
m4 = NC

!Part 4:
!------------------------------------------ Modifying Neibours of Neibours
Call setNeibour(Dim,Corn,Neib,N1,P1,NM,m1)
Call setNeibour(Dim,Corn,Neib,N2,P2,NM,m2)
Call setNeibour(Dim,Corn,Neib,N3,P1,NK,m3)
Call setNeibour(Dim,Corn,Neib,N4,P2,NK,m4)

!Part 5:

NP = NP + 1  ! Defining new point Nn
Nn = NP
X(Nn) = x_value
Y(Nn) = y_value

!Part 6:

Corn(m1,I1) = NM
Corn(m1,I2) = p1
Corn(m1,I3) = Nn
Corn(m1,4) = 0
Neib(m1,I1) = m3
Neib(m1,I2) = m2
Neib(m1,I3) = N1
Neib(m1,4) = 0

Corn(m2,I1) = NM
Corn(m2,I2) = Nn
Corn(m2,I3) = p2
Corn(m2,4) = 0
Neib(m2,I1) = m4
Neib(m2,I2) = N2
Neib(m2,I3) = m1
Neib(m2,4) = 0

Corn(m3,J1) = NK
Corn(m3,J2) = Nn
Corn(m3,J3) = p1
Corn(m3,4) = 0
Neib(m3,J1) = m1
Neib(m3,J2) = N3
Neib(m3,J3) = m4
Neib(m3,4) = 0

Corn(m4,J1) = NK
Corn(m4,J2) = p2
Corn(m4,J3) = Nn
Corn(m4,4) = 0
Neib(m4,J1) = m2
Neib(m4,J2) = m3
Neib(m4,J3) = N4
Neib(m4,4) = 0

print*,'Corn-m1',m1,':',Corn(m1,1),Corn(m1,2),Corn(m1,3)
print*,'Neib-m1',m1,':',Neib(m1,1),Neib(m1,2),Neib(m1,3)

print*,'Corn-m2',m2,':',Corn(m2,1),Corn(m2,2),Corn(m2,3)
print*,'Neib-m2',m2,':',Neib(m2,1),Neib(m2,2),Neib(m2,3)

print*,'Corn-m3',m3,':',Corn(m3,1),Corn(m3,2),Corn(m3,3)
print*,'Neib-m3',m3,':',Neib(m3,1),Neib(m3,2),Neib(m3,3)

print*,'Corn-m4',m4,':',Corn(m4,1),Corn(m4,2),Corn(m4,3)
print*,'Neib-m4',m4,':',Neib(m4,1),Neib(m4,2),Neib(m4,3)

!Part 7:

Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,m1)
Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,m2)
Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,m3)
Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,m4)

!===========================================================================================
End Subroutine AddNewPoint
!*********************************************************************************************