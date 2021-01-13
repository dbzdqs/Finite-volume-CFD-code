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
!// Developed by: H. Morad Tabrizi, Mechanical Eng., Amirkabir University of Technology    //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CST_InversToShap(NPtCurv_Up,NPtCurv_Lw,BSOrder_Up,BSOrder_Lw,ShapeFunc_Coe,&
                             Zita_te,X_Up,X_Lw,Ynew_Up,Ynew_Lw)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::NPtCurv_Up,NPtCurv_Lw,BSOrder_Up,BSOrder_Lw,ShapeFunc_Coe,Zita_te,X_Up,X_Lw
 Intent(Out  )::Ynew_Up,Ynew_Lw

 Integer::I,J,NPtCurv_Up,NPtCurv_Lw,BSOrder_Up,BSOrder_Lw,FN
 Real(8)::Zita_te,Choose,ShapeFun_Up,ShapeFun_Lw
 Real(8),Dimension(1:NPtCurv_Up)::X_Up,Ynew_Up
 Real(8),Dimension(1:NPtCurv_Lw)::X_Lw,Ynew_Lw
 Real(8),Dimension(1:BSOrder_Up+1+BSOrder_Lw+1)::ShapeFunc_Coe
!*********************************************************************************************
!Part1:
 Do I=1,NPtCurv_Up
 
  !Part2:	
   ShapeFun_Up = 0

  !Part3;
   Do J = 0,BSOrder_Up
      
	  ShapeFun_Up= ShapeFun_Up + ShapeFunc_Coe(J+1)* (X_Up(I)**(J)) *((1-X_Up(I))**(BSOrder_Up-J))* Choose( BSOrder_Up, J )

   End Do

  !Part4:
   Ynew_Up(I) = SQRT(X_Up(I))*(1-X_Up(I))* ShapeFun_Up + X_Up(I)* Zita_te

 End Do

!Part5:
 Do I=1,NPtCurv_Lw
  	
   ShapeFun_Lw = 0

   Do J = 0,BSOrder_Lw
      
	  ShapeFun_Lw= ShapeFun_Lw + ShapeFunc_Coe(BSOrder_Up+1+J+1)* (X_Lw(I)**(J)) *((1-X_Lw(I))**(BSOrder_Lw-J))* Choose( BSOrder_Lw, J )

   End Do

   Ynew_Lw(I) = SQRT(X_Lw(I))*(1-X_Lw(I))* ShapeFun_Lw + X_Lw(I)* Zita_te

 End Do
!*********************************************************************************************
 End
!###########################################################################################