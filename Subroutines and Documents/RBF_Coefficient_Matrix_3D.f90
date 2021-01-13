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
!// Date: Feb., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_Coefficient_Matrix_3D(Dim,X,Y,Z,List,nList,A)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,List,nList,X,Y,Z
 Intent(InOut)::A

 Integer::Dim,I,J,P1,P2,nList
 Real(8)::Dx,Dy,Dz,DL,Temp
 Integer,Dimension(1:Dim)::List
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:nList,1:nList)::A
!********************************************************************************************* 
!Part 1:
 Do I=1,nList
    Do J=I+1,nList

      !Part 2:
	   P1 = List(I)
	   P2 = List(J)

       Dx = X(P2)-X(P1)
       Dy = Y(P2)-Y(P1)
       Dz = Z(P2)-Z(P1)
       DL = Dsqrt( Dx*Dx + Dy*Dy + Dz*Dz )

      !Part 3:
       Call RBF_Function3D(DL,Temp)
       A(I,J)=Temp 

    End Do
 End Do

!Part 4:
 Do I=1,nList

   !Part 5:
    DL=0.0
    
   !Part 6:
	Call RBF_Function3D(DL,Temp)
    A(I,I)=Temp       

 End Do

!Part 7:
 Do I=1,nList
   Do J=I+1,nList
      A(J,I)=A(I,J) 
   End Do
 End Do
!*********************************************************************************************
 End
!###########################################################################################