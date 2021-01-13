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
!// Developed by: M. A. Zoljanahi, Mechanical Eng., Amirkabir University of Technology     //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_GreedyMovingMesh3D(Istp,Dim,NBP,NP,IBP,X,Y,Z,DelX,DelY,DelZ,NSBP,ACTV)
 Implicit None
!********************************************************************************************* 
 Intent(In    ):: Dim,NBP,NP,IBP,X,Y,Z,Istp
 Intent(inOut )::DelX,DelY,DelZ
 Intent(Out )::NSBP,ACTV

 Integer::Dim,I,NBP,NP,j,NSBP,Istp
 Real(8)::betax,betay,betaz,t2,t1
 Real(8),Dimension(1:Dim)::X,Y,Z,DelX,DelY,DelZ
 Integer,Dimension(1:Dim)::IBP
 Integer,Dimension(1:NBP)::ACTV
 Real(8),Dimension(1:NBP)::Cox,Coy,Coz,b_x,b_y,b_z
 Real(8),Dimension(1:NBP)::c_x,c_y,c_z,coxx,coyy,cozz
 Real(8),Dimension(1:NBP,1:NBP)::A
 Real(8),Dimension(1:NBP,1:NBP)::B
!********************************************************************************************* 
!Part 1:
 If (Istp==1) Then
     
    Call hybrid_greedy_algorithm_3D_2(NBP,Dim,X,Y,Z,IBP,Delx,Dely,Delz,NSBP,ACTV,coxx,coyy,cozz)  
  
 Else
     
    !part 2 
     Do I=1,NSBP
         
         c_x(I)=DelX(ACTV(I))
         c_y(I)=DelY(ACTV(I))
         c_z(I)=DelZ(ACTV(I))
         
     End Do
     
    !part 3  
     Call RBF_Coefficient_Matrix_3D(Dim,X,Y,Z,ACTV,NSBP,B)
     
    !part 4
     Call solve_lu(NSBP,B,c_x,coxx)
     Call solve_lu(NSBP,B,c_y,coyy)
     Call solve_lu(NSBP,B,c_z,cozz)
     
 End If
 
!part 5:
 Call RBF_Point_Dispalcement_3D(NSBP,Dim,ACTV,NP,X,Y,Z,DelX,DelY,DelZ,Coxx,Coyy,Cozz)
 
!*********************************************************************************************
 End
!###########################################################################################