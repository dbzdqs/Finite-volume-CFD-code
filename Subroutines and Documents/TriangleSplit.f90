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
Subroutine TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,tri,Nn)
Implicit None
!===========================================================================================
Intent(In)::Dim,Fronts,tri,States
Intent(Out)::Nn
Intent(InOut)::Corn,Neib,X,Y,NC,NP,FrontEdges

Integer::Dim,Fronts,A,B,C,N1,N2,N3,tri,NP,NC,Nn,t1,t2,t3,I
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
A = Corn(tri,1)
B = Corn(tri,2)
C = Corn(tri,3)

N1 = Neib(tri,1)
N2 = Neib(tri,2)
N3 = Neib(tri,3)

NP = NP + 1
Nn = NP

X(Nn) = (X(A) + X(B) + X(C))/3
Y(Nn) = (Y(A) + Y(B) + Y(C))/3
!Part 2:
t1 = tri
NC = NC + 1
t2 = NC
NC = NC + 1
t3 = NC

Corn(t1,1) = A
Corn(t1,2) = B
Corn(t1,3) = Nn
Neib(t1,1) = t3
Neib(t1,2) = t2
Neib(t1,3) = N3

Corn(t2,1) = A
Corn(t2,2) = Nn
Corn(t2,3) = C
Neib(t2,1) = t3
Neib(t2,2) = N2
Neib(t2,3) = t1

Corn(t3,1) = Nn
Corn(t3,2) = B
Corn(t3,3) = C
Neib(t3,1) = N1
Neib(t3,2) = t2
Neib(t3,3) = t1
!Part 3:
do I=1,4
	if(Neib(N1,I) == tri) Neib(N1,I) = t3
	if(Neib(N2,I) == tri) Neib(N2,I) = t2
	if(Neib(N3,I) == tri) Neib(N3,I) = t1
end do
!Part 4:
Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,t1)
Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,t2)
Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,t3)
!===========================================================================================
End Subroutine TriangleSplit
!*********************************************************************************************
