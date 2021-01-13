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
!// Developed by: F. Farhadkhani, Mathmatical, Amirkabir university of Technology          //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Find_PointEdge(Dim,EI,IDS,InxEdgeOfCell,EP1,EP2,E1,E2,E3,E4)
Implicit None
!*********************************************************************************************
Intent(In   )::Dim,EI,IDS,InxEdgeOfCell
Intent(Out  )::EP1,EP2,E1,E2,E3,E4

Integer::Dim,I,ME,NE,P1,P2,N,EI,EP1,EP2,E1,E2,E3,E4
Integer,Dimension(1:4,1:Dim)::IDS
Integer,Dimension(1:4,1:Dim)::InxEdgeOfCell
!*********************************************************************************************
!Part 1:
 ME = IDS(1, EI)
 NE = IDS(2, EI)
 P1 = IDS(3, EI)
 P2 = IDS(4, EI)

!Part 2:
 Do I=1,3
 
   !Part 3:
    N = abs(InxEdgeOfCell(I,ME))

   !Part 4:
	If(N/=EI)Then

    !Part 5:
     If(IDS(3, N) == P1) Then
	  E1 =  N
      EP1 =IDS(4, N)
     Elseif(IDS(4, N) == P1) Then
	  E1 =  N
      EP1 =IDS(3, N)
	 Endif

    !Part 6:
     If(IDS(3, N) == P2 .or. IDS(4, N) == P2) E3 =  N

    Endif

 End do

!Part 7:
 Do I=1,3

    N = abs(InxEdgeOfCell(I,NE))

	If(N/=EI)Then

     If(IDS(3, N) == P1) Then
	  E4 =  N
      EP2 =IDS(4, N)
     Elseif(IDS(4, N) == P1) Then
	  E4 =  N
      EP2 =IDS(3, N)
	 Endif

     If(IDS(3, N) == P2 .or. IDS(4, N) == P2) E2 =  N

    Endif

 End do
!*********************************************************************************************
 End
!###########################################################################################
    
