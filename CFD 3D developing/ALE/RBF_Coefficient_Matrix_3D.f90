!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Coefficient Matrix for Radfial Basis Function             //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: Agust, 03, 2015                                                                //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Developed by: A.R. Rezaei, Iran, Tehran, a.r.rezaei@aut.ac.ir                        //!
!// Doc ID: MC5F029F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine RBF_Coefficient_Matrix_3D(Dim,X,Y,Z,List,nList,A)
 Implicit None
!******************************************************************************************* 
 Intent(In   )::Dim,List,nList,X,Y,Z
 Intent(InOut)::A

 Integer::Dim,I,J,P1,P2,nList
 Real(8)::Dx,Dy,Dz,DL,Temp
 Integer,Dimension(1:Dim)::List
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:nList,1:nList)::A
!******************************************************************************************* 
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
       Call RBF_Function(DL,Temp)
       A(I,J)=Temp 

    End Do
 End Do

!Part 4:
 Do I=1,nList

   !Part 5:
    DL=0.0
    
   !Part 6:
	Call RBF_Function(DL,Temp)
    A(I,I)=Temp       

 End Do

!Part 7:
 Do I=1,nList
   Do J=I+1,nList
      A(J,I)=A(I,J) 
   End Do
 End Do
!*******************************************************************************************
 End
!###########################################################################################