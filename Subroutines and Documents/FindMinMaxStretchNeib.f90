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
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine FindMinMaxStretchNeib(Dim,NP,Vertice,X,Y,Z,NConnectedPoints,IConnectedPoints,NodeMetric,NExceptNodes,IExceptNodes,MinNodeIndex,MaxNodeIndex,MinValue,MaxValue)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,MinNodeIndex,MaxNodeIndex,I,J,Vertice,NExceptNodes,NConnectedNodes,NConnectedPoints
 Logical::existsInExceptList
 Real(8),Dimension(1:3,1:3,1:Dim)::NodeMetric
 Real(8)::hV,MinValue,MaxValue
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(100)::IExceptNodes
 Integer,Dimension(1:100)::IConnectedPoints
!*********************************************************************************************
!Part 1:
 MinValue    = 0.0
 MaxValue    = 0.0
 MinNodeIndex = 0
 MaxNodeIndex = 0
 Do I=1,NConnectedPoints
     
   !Part 2:
    existsInExceptList = .FALSE.
    Do J=1,NExceptNodes
       If(IExceptNodes(J)==IConnectedPoints(I))Then
        existsInExceptList = .TRUE.
        exit
       Endif
    EndDo
     
   !Part 3:
    Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,Vertice),Vertice,IConnectedPoints(I),hV)
    If( (hV<MinValue .OR. MinNodeIndex==0) .AND. (.NOT. existsInExceptList) )Then
     MinValue    = hV
     MinNodeIndex = IConnectedPoints(I)
    Endif
    If( (hV>MaxValue .OR. MaxNodeIndex==0) .AND. (.NOT. existsInExceptList) )Then
     MaxValue    = hV
     MaxNodeIndex = IConnectedPoints(I)
    Endif
     
 EndDo
!*********************************************************************************************
End