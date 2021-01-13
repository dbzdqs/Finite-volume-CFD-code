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
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
  Subroutine FindMinStretchNeib(Dim,NP,Vertice,NConnectedNodes,IConnectedNodes,NodeMetric,NExceptNodes,IExceptNodes,NodeIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,Node1,Node2,NodeIndex,I,J,Vertice
 Logical::existsInExceptList
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8)::h12,h21,Norm,Length,EucLen,minVal,hV
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1,2)::UnitV
 Real(8),Dimension(2,1)::UnitVt
 Integer,Dimension(1:2)::MDim1,MDim2
 Real(8),Dimension(1)::TempMatrix1,TempMatrix4
 Real(8),Dimension(1:2,1)::TempMatrix2
 Real(8),Dimension(1:2,1:10)::Vector
 Real(8),Dimension(1:Dim)::ElipsoidRatio
 Integer::NExceptNodes
 Integer,Dimension(200)::IExceptNodes
 Integer,Dimension(1:100)::IConnectedNodes
 Integer::NConnectedNodes
!*********************************************************************************************
!Part 1:
 minVal    = 0.0
 NodeIndex = 0
 Do I=1,NConnectedNodes
     
   !Part 2:
    existsInExceptList = .FALSE.
    Do J=1,NExceptNodes
       If(IExceptNodes(J)==IConnectedNodes(I))Then
        existsInExceptList = .TRUE.
        exit
       Endif
    EndDo
     
   !Part 3:
    Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,Vertice),Vertice,IConnectedNodes(I),hV)
    If( (hV<minVal .OR. NodeIndex==0) .AND. (.NOT. existsInExceptList) )Then
     minVal    = hV
     NodeIndex = IConnectedNodes(I)
    Endif
     
 EndDo
!*********************************************************************************************
End