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
!// Date: Aug., 30, 2015                                                                   //!
!// Developed by: A. Moslemi Pak, Mechanical Eng., Amirkabir University of Technology      //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CalCellCentroid(Dim,NF,NC,IDS,X,Y,XCen,YCen)
 Implicit None
!*********************************************************************************************
 Intent(In   )                   ::  Dim,NF,NC,IDS,X,Y
 Intent(Out  )                   ::  XCen,YCen

 Integer                         ::  Dim,NF,NC,i,J,ME,NE,P1,P2
 Real(8)                         ::  X1,Y1,X2,Y2
 Integer,Dimension(1:4,1:Dim)    ::  IDS
 Real(8),Dimension(1:Dim)        ::  X,Y,XCen,YCen
 Integer,Dimension(1:Dim)        ::  NPCell
!*********************************************************************************************
!Part 1: Predefining the centroid of the cells coordinates
 XCen(1:NC)   = 0.0
 YCen(1:NC)   = 0.0
 NPCell(1:NC) = 0
    
!Part 2: Finding the Centroid
 Do i=1,NF
    
    ME = IDS(1,i)
    NE = IDS(2,i)
    P1 = IDS(3,i)
    P2 = IDS(4,i)
    
    X1 = X(P1)
    Y1 = Y(P1)
    
    X2 = X(P2)
    Y2 = Y(P2)
    
   !The First Cell of the Face
    XCen(ME)      = XCen(ME) + X1 + X2
    YCen(ME)      = YCen(ME) + Y1 + Y2
    NPCell(ME)    = NPCell(ME) + 2

   !The Second Cell of the Face
    If (NE/=0) Then
        XCen(NE)   = XCen(NE) + X1 + X2
        YCen(NE)   = YCen(NE) + Y1 + Y2
        NPCell(NE) = NPCell(NE) + 2
    End If
    
 End Do

!Part 3: Finalize the center of the elements
 Do i=1,NC
    XCen(i) = XCen(i)/NPCell(i)
    YCen(i) = YCen(i)/NPCell(i)
 End Do
!*********************************************************************************************
 End Subroutine CalCellCentroid
!###########################################################################################