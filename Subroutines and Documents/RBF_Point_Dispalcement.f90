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
 Subroutine RBF_Point_Dispalcement(Dim,IBP,NP,NBP,Cox,Coy,X,Y,DelX,DelY)
 Implicit None
!********************************************************************************************* 
 Intent(In ):: Dim,IBP,NP,NBP,Cox,Coy,X,Y
 Intent(InOut )::DelX,DelY

 Integer::Dim,I,J,NP,NBP
 Real(8)::Dx,Dy,DL,Temp,Sum_DelX,Sum_DelY,Xi,Yi
 Integer,Dimension(Dim)::IBP
 Real(8),Dimension(Dim)::X,Y,DelX,DelY
 Real(8),Dimension(NBP+3)::Cox,Coy
!********************************************************************************************* 
!Part 1:
 Do I=1,NP

   !Part 2:
    Xi = X(I)
	Yi = Y(I)

   !Part 3:
    Sum_DelX = 0.0
	Sum_DelY = 0.0

   !Part 4:            
    Do J=1,NBP

      !Part 5:
	   IF(I==IBP(J)) goto 10

      !Part 6:
       Dx = X( IBP(J) )-Xi
       Dy = Y( IBP(J) )-Yi
       DL = Dsqrt( Dx*Dx + Dy*Dy )
      
      !Part 7:        
       Call RBF_Function(DL,Temp)
    
      !Part 8:     
       Sum_DelX = Sum_DelX + Cox(J)*Temp   
       Sum_DelY = Sum_DelY + Coy(J)*Temp

    End Do
   
   !Part 9:
    DelX(I) = Sum_DelX + Cox(NBP+1)*1. + Cox(NBP+2)*Xi + Cox(NBP+3)*Yi   
    DelY(I) = Sum_DelY + Coy(NBP+1)*1. + Coy(NBP+2)*Xi + Coy(NBP+3)*Yi

10 End Do
!*********************************************************************************************
 End
!###########################################################################################
    
