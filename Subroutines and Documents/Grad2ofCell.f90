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
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Grad2AtCell(Dim,NC,NF,NFW1,NF1,NF2,IDS,A,NX,NY,DUY_C,DD2UY)         
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF,NFW1,NF1,NF2,IDS,A,NX,NY,DUY_C
 Intent(Out  )::DD2UY

 Integer::Dim,I,NC,NFW1,NF,NF1,NF2,ME,NE
 Real(8)::DU
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::A
 Real(8),Dimension(1:Dim)::DUY_C,DD2UY,NX,NY
!***************************************************************************************************  
!Part 1:
 DO I=1,NC    
    DD2UY(I) = 0.0
 End Do

!Part 2:
 DO I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 4:
    DU  = 0.5*( DUY_C(ME)+DUY_C(NE) )

   !Part 5:
    DD2UY(ME) = DD2UY(ME) + DU*NY(I)
    DD2UY(NE) = DD2UY(NE) - DU*NY(I)
 End Do

!Part 6:
 DO I=NFW1+1,NF

    ME = IDS(1,I)
    DU  = DUY_C(ME)
    DD2UY(ME) = DD2UY(ME) + DU*NY(I)
    
 End Do

!Part 7:
 DO I=1,NC
    DD2UY(I) = DD2UY(I)/A(I)
 End Do
!*********************************************************************************************
 End
!###########################################################################################