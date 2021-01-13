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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_Point_Dispalcement_3D(nList,Dim,List,NP,X,Y,Z,DelX,DelY,DelZ,cx,cy,cz)
 Implicit None
!********************************************************************************************* 
 Intent(In ):: Dim,List,NP,X,Y,Z,cx,cy,cz,nList
 Intent(Out)::DelX,DelY,DelZ
  
 Integer::Dim,I,J,NP,nList
 Real(8)::Dx,Dy,Dz,DL,Temp,Sum_DelX,Sum_DelY,Sum_DelZ,Xi,Yi,Zi
 Integer,Dimension(Dim)::List
 Real(8),Dimension(Dim)::X,Y,Z,DelX,DelY,DelZ
 Real(8),Dimension(nList)::cx,cy,cz
!********************************************************************************************* 
!Part 1:
 Do I=1,NP

   !Part 2:
    Xi = X(I)
	Yi = Y(I)
    Zi = Z(I)

   !Part 3:
    Sum_DelX = 0.0
	Sum_DelY = 0.0
    Sum_DelZ = 0.0

   !Part 4:            
    Do J=1,nList

      !Part 5:
	   IF(I==List(j)) goto 8

      !Part 6:
       Dx = X( List(j) )-Xi
       Dy = Y( List(j) )-Yi
       Dz = Z( List(j) )-Zi
       
       DL = Dsqrt( Dx*Dx + Dy*Dy + Dz*Dz)
      
      !Part 7:        
       Call RBF_Function3D(DL,Temp)
    
      !Part 8:     
       Sum_DelX = Sum_DelX + cx(j)*Temp   
       Sum_DelY = Sum_DelY + cy(j)*Temp
       Sum_DelZ = Sum_DelZ + cz(j)*Temp
    End Do
   
   !Part 9:
    DelX(I) = Sum_DelX   
    DelY(I) = Sum_DelY
    DelZ(I) = Sum_DelZ

8 End Do
!*********************************************************************************************
 End
!###########################################################################################
    
