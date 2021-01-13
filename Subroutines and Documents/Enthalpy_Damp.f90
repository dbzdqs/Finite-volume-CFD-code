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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Enthalpy_Damp(Dim,NC,GM,DT,H0,Wnp1)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,GM,DT,H0
 Intent(InOut)::Wnp1

 Integer::J,NC,Dim,Ncyc
 Real(8)::U,V,GM,H0,Alpha,del,Beta,C,Mach
 Real(8),Dimension(1:Dim)::P,H,DT
 Real(8),Dimension(1:4,1:Dim)::Wnp1
!*********************************************************************************************	
!Part 1:
 Beta = 0.0
 
 Do J=1,NC

   !Part 2:
    Alpha = 0.15

   !Part 3:
    U    = Wnp1(2,J)/Wnp1(1,J)
    V    = Wnp1(3,J)/Wnp1(1,J)
    P(J) = (GM-1)*(Wnp1(4,J)-0.5*Wnp1(1,J)*(U*U+V*V))

   !Part 4:
    Mach = DSqrt( (U*U+V*V) / (GM*P(J)/Wnp1(1,J)) )

   !Part 5:
    If(Mach>1.) Alpha=0.0

   !Part 6:
    H(J) = GM/(GM-1)*P(J)/Wnp1(1,J)+0.5*(U*U+V*V)

   !Part 7:
    C   = 1+Alpha*DT(J)*(H(J)-H0)
    Del = Beta*DT(J)*Wnp1(1,J)*(H(J)-H0)/GM

   !Part 8: 
    Wnp1(1,J) = Wnp1(1,J) / C
    Wnp1(2,J) = Wnp1(2,J) / C
    Wnp1(3,J) = Wnp1(3,J) / C
    Wnp1(4,J) = (Wnp1(4,J)-Del) / C
 end Do
!*********************************************************************************************
 End
!###########################################################################################