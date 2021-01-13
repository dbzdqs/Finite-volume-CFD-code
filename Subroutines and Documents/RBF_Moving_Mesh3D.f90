!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate Displacement of non-Boundary Points by Radfial Basis Function //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: Agust, 03, 2015                                                                //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Developed by: A.R. Rezaei, Iran, Tehran, a.r.rezaei@aut.ac.ir                        //!
!// Doc ID: MC5F031F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine RBF_Moving_Mesh3D(Dim,NBP,NP,IBP,X,Y,Z,DelX,DelY,DelZ)
 Implicit None
!******************************************************************************************* 
 Intent(In    ):: Dim,NBP,NP,IBP,X,Y,Z
 Intent(InOut )::DelX,DelY,DelZ

 Integer::Dim,I,N,NBP,NP,j
 Real(8),Dimension(1:Dim)::X,Y,Z,DelX,DelY,DelZ
 Integer,Dimension(1:Dim)::IBP
 Real(8),Dimension(1:NBP)::Cox,Coy,Coz,b_x,b_y,b_z
 Real(8),Dimension(1:NBP,1:NBP)::A
!******************************************************************************************* 
!Part 1:
 Do I=1,NBP
    N = IBP(I)

    b_x(I)=DelX(N)
    b_y(I)=DelY(N)
    b_z(I)=DelZ(N)
 End Do
 
!Part 2:
 Call RBF_Coefficient_Matrix_3D(Dim,X,Y,Z,IBP,NBP,A)

!Part 3:
 call solve_lu(NBP,A,b_x,cox)
 call solve_lu(NBP,A,b_y,coy)
 call solve_lu(NBP,A,b_z,coz)

!Part 4:
 Call RBF_Point_Dispalcement_3D(NBP,Dim,IBP,NP,X,Y,Z,DelX,DelY,DelZ,Cox,Coy,Coz)
!*******************************************************************************************
 End
!###########################################################################################