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
!// Developed by: M. Vakili, Computer Science, Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine WriteMsh(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim)
 Implicit None
!*********************************************************************************************
 Intent(In)::Dim,DimIDS ,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim

 Integer::Dim,DimIDS ,NP,NC,NR,NF,I,J,L,FaceCount,MeshDim
 Integer,Dimension(1:10)::NFR,BCType
 Integer,Dimension(1:Dim)::FaceType,CellType
 Integer,Dimension(1:DimIDS,1:Dim)::IDS
 Real*8,Dimension(1:Dim)::X,Y,Z
 Character*100,Dimension(1:10)::BCTitle
!*********************************************************************************************
!Part1
 Open(1,File='MeshOut.msh')

!Part2
 Write(1,'(A)') '(0 "This Mesh is generated by an open source code")'

!Part3	
 Write(1,'(A,I0,A)') '(2 ',MeshDim,')'

!Part4    
 Write(1,'(A)') '(0 "Node Section")'
    
!Part5
 Write(1,'(A,Z0,A,I0,A)') '(10 (0 1 ',NP,' 0 ',MeshDim,'))'
    
!Part6
 Write(1,'(A,Z0,A,I0,A)') '(10 (3 1 ',NP,' 1 ',MeshDim,')'
    
!Part7
 Write(1,'(A)') '('

!Part8
 If(MeshDim==2) Then
  Do I=1,NP
     Write(1,'(1X,ES17.10,3X,ES17.10)') X(I),Y(I)
  End Do

!PArt9
 Else
  Do I=1,NP
     Write(1,'(1X,ES17.10,3X,ES17.10,3X,ES17.10)') X(I),Y(I),Z(I)
  End Do
 End If
 
!Part10    
 Write(1,'(A)') '))'

!Part11    
 Write(1,'(A,Z0,A)') '(12 (0 1 ',NC,' 0 0))'

!Part12
 Write(1,'(A,Z0,A)') '(12 (4 1 ',NC,' 1 0)'

!Part13
 Write(1,'(A1)', Advance='No') '('
 Do J=1,NC
    Write(1,'(I2)', Advance='No') CellType(J)
    If(Mod(J,20) == 0) Then
     Write(1,*)
    End If
 End Do
 Write(1,*)
 Write(1,'(A)') '))'

!Part14    
 Write(1,'(A,Z0,A)') '(13 (0 1 ',NF,' 0 0))'
     
!Part15
 FaceCount=0
 Do I=1,NR

   !Part16
    If(Len_Trim(BCTitle(I))>0) Then
     Write(1,'(A,A,A)') '(0 "',Trim(AdjustL(BCTitle(I))),'")'
    End If

   !Part17
    Write(1,'(A,Z0,1X,Z0,1X,Z0,1X,Z0,A,Z0,A)') '(13 (',4+I,FaceCount+1,FaceCount+NFR(I),BCType(I),' 0)('
        
   !Part18
    Do J=FaceCount+1,FaceCount+NFR(I)

      !Part19
       Write(1,'(I1,1X)', Advance='No') FaceType(J)

      !Part20
       Do L=3,FaceType(J)+2            
           Write(1,'(Z0,1X)', Advance='No') IDS(L,J)
       End Do

      !Part21
       Write(1,'(Z0,1X)', Advance='No') IDS(1,J)
       Write(1,'(Z0,1X)', Advance='No') IDS(2,J)
       Write(1,*)           
            
    End Do

   !Part22
    Write(1,'(A)') ')'
    Write(1,'(A)') ')'

   !Part23
    FaceCount=FaceCount+NFR(I)
        
 End Do	

 Close(1)
!*********************************************************************************************
 End Subroutine
!###########################################################################################
