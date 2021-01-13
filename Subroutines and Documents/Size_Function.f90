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
 Subroutine Size_Function(Dim,BFacPt,NFacCrv,NBoundCrv,X,Y,SF)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,BFacPt,NFacCrv,NBoundCrv,X,Y
 Intent(Out  )::SF

 Integer::Dim,J,J1,P1,P2,SumE,NBoundCrv
 Real(8)::Edge_L
 Real(8),Dimension(1:Dim)::X,Y,SF
 Integer,Dimension(1:Dim,1:2)::BFacPt
 Integer,Dimension(1:10)::NFacCrv
!*********************************************************************************************
!Part 1:
 SF=0.0

!Part 2:
 SumE=0
 Do J=1,NBoundCrv 
     Do J1=SumE+1,SumE+NFacCrv(J)

       !Part 3:
        P1 = BFacPt(J1,1)
        P2 = BFacPt(J1,2)

       !Part 4:
	    Edge_L = Sqrt( (X(P1)-X(P2))**2 + (Y(P1)-Y(P2))**2 )

       !Part 5:
	    SF(P1) = SF(P1) + Edge_L / 2
	    SF(P2) = SF(P2) + Edge_L / 2

    End Do
    SumE = SumE + NFacCrv(J)
 End Do
!*********************************************************************************************
 End 
!###########################################################################################
