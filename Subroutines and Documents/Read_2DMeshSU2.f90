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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: A. Moslemi Pak, Mechanical Eng., Amirkabir University of Technology      //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Read_2DMeshSU2(Dim,MeshDim,NC,Corn,NP,X,Y,NBoundCrv,BCName,NFacCrv,BFacPt,CellType)
    Implicit None
!*********************************************************************************************
    Intent(In   )                           ::  Dim
    Intent(Out  )                           ::  MeshDim,NC,Corn,NP,X,Y,NBoundCrv,BCName,NFacCrv,BFacPt,CellType

    Integer                                 ::  Dim,IO,MeshDim,NC,i,j,NP,NBoundCrv,NBE,Dumy,Dumy2
    Character*100                           ::  Temp
    Character*100,Dimension(1:10)           ::  BCName
    Integer,Dimension(1:4,1:Dim)            ::  Corn
    Integer,Dimension(1:100)                ::  NFacCrv
    Real(8),Dimension(1:Dim)                ::  X,Y
    Integer,Dimension(1:Dim)                ::  CellType
    Integer,Dimension(1:Dim,1:2)            ::  BFacPt
    Integer,Dimension(1:3)                  ::  TempBFacPt
!*********************************************************************************************
    !Part 1: Opening the SU2 Mesh File
    Open(1,File='Mesh.su2')
    
    !Part 2: Read dimension of the space
    Read(1,'(A100)',IOSTAT = IO) Temp
    Do While (Temp(1:5) /= 'NDIME')
        Read(1,'(A100)',IOSTAT = IO) Temp
    End Do
    Read(Temp(7:100),*) MeshDim
    
    !Part 3: Read Number of the elements
    Read(1,'(A100)',IOSTAT = IO) Temp
    Do While (Temp(1:5) /= 'NELEM')
        Read(1,'(A100)',IOSTAT = IO) Temp
    End Do
    Read(Temp(7:100),*) NC
    
    !Part 4: The node number connectivity
    Corn(1:4,1:NC) = 0
    Do i=1,NC
        Read(1,'(A2)',ADVANCE='NO') Temp
        Read(Temp(1:2),*) Dumy2
        !Check the 'VTK' format
        If (Dumy2 == 5) Then
            !This is a triangular element and it contains only three nodes
            Read(1,*) Corn(1,i), Corn(2,i), Corn(3,i)
            !Adding 1 for fitting the format for the tecplot
            Corn(1:3,i) = Corn(1:3,i) + 1
            CellType(i) = 3
        Else If (Dumy2 == 9) Then
            !This is a quadrilateral element and it contains only four nodes
            Read(1,*) Corn(1,i), Corn(2,i), Corn(3,i), Corn(4,i)
            !Adding 1 for fitting the format for the tecplot
            Corn(1:4,i) = Corn(1:4,i) + 1
            CellType(i) = 4
        Else If (Dumy2 == 3) Then
            !This is a line element and it contains only two nodes
            Read(1,*) Corn(1,i), Corn(2,i)
            !Adding 1 for fitting the format for the tecplot
            Corn(1:2,i) = Corn(1:2,i) + 1
            CellType(i) = 0
        End If
    End Do
    
    !Part 5: Read the Points Number
    Read(1,'(A100)',IOSTAT = IO) Temp
    Do While (Temp(1:5) /= 'NPOIN')
        Read(1,'(A100)',IOSTAT = IO) Temp
    End Do
    Read(Temp(7:100),*) NP
    
    !Part 6: Read the Points Coordinates
    Do i=1,NP
        Read(1,*) X(i) , Y(i)
    End Do
    
    !Part 7: Reading the Boundaries data
    Read(1,'(A100)',IOSTAT = IO) Temp
    Do While (Temp(1:5) /= 'NMARK')
        Read(1,'(A100)',IOSTAT = IO) Temp
    End Do
    Read(Temp(7:100),*) NBoundCrv
    
    !Part 8: Reading the boundary curves points
    NBE = 0
    Do i=1,NBoundCrv
        !Read the boundary name
        Read(1,'(A100)',IOSTAT = IO) Temp
        Read(Temp(13:100),*) BCName(i)
        !Read the number of parts of each boundary
        Read(1,'(A100)',IOSTAT = IO) Temp
        Read(Temp(15:100),*) NFacCrv(i)
        TempBFacPt(1:3)=0
        Do j=(1+NBE),(NBE+NFacCrv(i))
            Read(1,*) TempBFacPt(1), TempBFacPt(2), TempBFacPt(3)
            BFacPt(j,1:2) = TempBFacPt(2:3) + 1
        End Do
        NBE = NBE + NFacCrv(i)
    End Do
    
    Close(1)
!*********************************************************************************************    
End Subroutine Read_2DMeshSU2