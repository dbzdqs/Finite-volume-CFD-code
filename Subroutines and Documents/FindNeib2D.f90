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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: A. Moslemi Pak, Mechanical Eng., Amirkabir University of Technology      //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine FindNeib2D(Dim,NC,Corn,CellType,Neib)
    Implicit None
!*********************************************************************************************
    Intent(In   )                           ::  Dim,NC,Corn,CellType
    Intent(Out   )                          ::  Neib

    Integer                                 ::  Dim,NC,j,j1,Count,q,Vertices1,Dumy2
    Integer                                 ::  Vertices2,m,n
    Integer,Dimension(1:4,1:Dim)            ::  Neib,Corn
    Integer,Dimension(1:5)                  ::  TempNeib1,TempNeib2
    Integer,Dimension(1:Dim)                ::  CellType
!*********************************************************************************************
    !Part 1:
    Neib(1:4,1:NC) = 0
    
    Do j=1,NC
        !Check the number of nodes containing the current cell
        Vertices1 = CellType(j)
        TempNeib1(1:(Vertices1+1))   = 0
        Count = 0
        Do q=1,Vertices1
            TempNeib1(q) = Corn(q,j)
        End Do
        TempNeib1(Vertices1+1) = Corn(1,j)
        Do j1=1,NC
            If (j==j1) GOTO 500
            !Check the number of nodes containing the current cell
            Vertices2 = CellType(j1)
            TempNeib2(1:(Vertices2+1))   = 0
            
            Do q=1,Vertices2
                TempNeib2(q) = Corn(q,j1)
            End Do
            TempNeib2(Vertices2+1) = Corn(1,j1)
            Do m=1,Vertices1
                Do n=1,Vertices2
                    If ( (TempNeib1(m)==TempNeib2(n) .AND. TempNeib1(m+1)==TempNeib2(n+1)) .OR.  &
                        ( (TempNeib1(m)==TempNeib2(n+1) .AND. TempNeib1(m+1)==TempNeib2(n)) ) ) Then
                        Dumy2 = m-1
                        If (Dumy2==0) Dumy2=Vertices1
                        Count = Count + 1
                        Neib(Dumy2,j) = j1
                        !Skip for the next element (cycling "j") when all of its neighbours were found completely
                        If (Count == Vertices1) GOTO 550
                        !Skip for the next element (cycling "j1") when it was the neighbour of element "j"
                        GOTO 500
                    End If
                End Do
            End Do
500     End Do
550     End Do
    
    
End Subroutine FindNeib2D
