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
!// Date: May., 15, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************  
 Subroutine Find_EdgePointNeib(Dim,Dead,IDS,InxEdgOfCell,EI,EP1,EP2,E1,E2,E3,E4,N1,N2,N3,N4)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,IDS,InxEdgOfCell,EI,Dead
 Intent(Out  )::EP1,EP2,E1,E2,E3,E4,N1,N2,N3,N4

 Integer::Dim,I,ME,NE,EP1,EP2,P1,P2,N,EI,E1,E2,E3,E4,N1,N2,N3,N4,EdgesCount,Dead
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell
!*********************************************************************************************
!Part 1:
 ME = IDS(1, EI)
 NE = IDS(2, EI)
 P1 = Dead
 If(IDS(3, EI)==P1)Then
     P2 = IDS(4, EI)
 Else
     P2 = IDS(3, EI)
 EndIf
 
!Part 2:
 EdgesCount = 3
 If(InxEdgOfCell(4,ME)/=0)  EdgesCount = 4

!Part 3:
 Do I=1,EdgesCount
     
   !Part 4:
    N = abs(InxEdgOfCell(I,ME))    
	
   !Part 5:
    If(N/=EI)Then
     If(IDS(3, N) == P1)Then
      E1    = N
      EP1   = IDS(4,N)
     Elseif(IDS(4, N) == P1)Then
      E1    = N
      EP1   = IDS(3,N)
     Endif
     
    !Part 6:
     If(IDS(3, N) == P2 .or. IDS(4, N) == P2)  E2 =  N
            
    Endif
 End do

!Part 7:
 If (IDS(1,E1)==ME) Then
  N1 = 2
 Else
  N1 = 1
 End If
 
!Part 8:
 If (IDS(1,E2)==ME) Then
  N2 = 1
 Else
  N2 = 2
 END IF

 If(NE/=0)Then
 !Part 9:
  EdgesCount = 3
  If(InxEdgOfCell(4,NE)/=0)  EdgesCount = 4
  
 !Part 10:
  Do I=1,EdgesCount
  
	!Part 11:
     N = abs(InxEdgOfCell(I,NE))
     If(N/=EI)Then
	 !Part 12:
      If(IDS(3, N) == P1) Then
       E3   = N
       EP2  = IDS(4, N)
      Elseif(IDS(4, N) == P1) Then
       E3   = N
       EP2  = IDS(3, N)
      Endif
      
	  !Part 13:
      IF(IDS(3, N) == P2 .or. IDS(4, N) == P2) E4 =  N
      
     Endif
  End do
 
 !Part 14:
  If (IDS(1,E3)==NE) Then
   N3 = 2
  Else
   N3 = 1
  END If

 !Part 15: 
  IF (IDS(1,E4)==NE) Then
   N4 = 1
  Else
   N4 = 2
  END IF
  
 END IF
!*********************************************************************************************
 End
!###########################################################################################