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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CalDisFunction(Dim,NBP,IBP,X,Y,CoRBF,Xp,Yp,DisFuncPt)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NBP,IBP,X,Y,CoRBF,Xp,Yp
 Intent(Out  )::DisFuncPt

 Integer::Dim,J,NBP
 Real(8) ::DisFuncPt,Dx,Dy,DL,Temp,P1
 Real(8),dimension(1:Dim)::X,Y
 Real(8)::CoRBF(1:NBP)
 Real(8):: XP,YP
 Integer,Dimension(Dim)::IBP 
!*********************************************************************************************
!Part 1:
 DisFuncPt = 0.0

!Part 2:            
 Do J=1,NBP

    P1 = IBP(J)
       
   !Part 3:
    Dx = X(P1)-XP
    Dy = Y(P1)-YP
    DL = Dsqrt( Dx*Dx + Dy*Dy )
      
   !Part 4:        
    Call RBF_FunctionV2(DL,Temp)

   !Part 5:     
    DisFuncPt = DisFuncPt + CoRBF(J)*Temp

 End Do
!*********************************************************************************************
 End 
!###########################################################################################