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
!// Date: Mar., 05, 2013                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConectedCellOfPoint(Dim,Point,IDS,NConectEdge,IConectEdge,NConectCell,IConectCell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Point,IDS,NConectEdge,IConectEdge
 Intent(Out  )::NConectCell,IConectCell

 Integer::Dim,I,J,ME,NE,E,Repeated,NConectCell,Point
 Integer,Dimension(1:Dim)::NConectEdge
 Integer,Dimension(1:Dim,1:100)::IConectEdge
 Integer,Dimension(1:50)::IConectCell
 Integer,Dimension(4,Dim)::IDS
!*********************************************************************************************
!Part 1:
 NConectCell = 0
 
!Part 2:
 Do J=1,NConectEdge(Point)

   !Part 3:
    E   = IConectEdge(Point,J)
    ME  = IDS(1,E)
    NE  = IDS(2,E)

   !Part 4:
    Repeated = -1
    Do I=1,NConectCell
       IF(IConectCell(I) == ME) Repeated = 1
    End do
	
   !Part 5:
    IF( Repeated==-1 )Then
     NConectCell = NConectCell + 1 
     IConectCell(NConectCell) = ME
	End IF

   !Part 6:
    Repeated = -1
    Do I=1,NConectCell
       IF( IConectCell(I)==NE .or. NE==0 ) Repeated = 1
    End Do
    
   !Part 7:
    IF( Repeated==-1 )Then
     NConectCell = NConectCell + 1 
     IConectCell(NConectCell) = NE
    End IF

 End Do
!*********************************************************************************************
 End
!###########################################################################################