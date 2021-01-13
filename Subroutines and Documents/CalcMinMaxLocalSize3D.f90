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
Subroutine CalcMinMaxLocalSize3D(Dim,NP,Metric,Node,X,Y,Z,NConnectedPoints,IConnectedPoints,MinVal,MaxVal,EDir)
 Implicit None
!*********************************************************************************************
 Integer::Dim,I,EDir,NP,Node,Node2
 Real(8),Dimension(1:3,1:3)::Metric
 Real(8)::MinVal,MaxVal,h
 Logical::flag1
 Integer::NConnectedPoints
 Integer,Dimension(1:100)::IConnectedPoints
 Real(8),Dimension(1:Dim)::X,Y,Z
!*********************************************************************************************
!Part 1:
 flag1=.TRUE.
 Do I=1,NConnectedPoints
    Node2=IConnectedPoints(I)
    Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,Metric,Node,Node2,h)
     
   !Part 2:
    If(MinVal>h .OR. flag1) MinVal = h
    If(MaxVal<h .OR. flag1)Then
        MaxVal = h
        EDir   = I
    Endif
    flag1 = .FALSE.
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################