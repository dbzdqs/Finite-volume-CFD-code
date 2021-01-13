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
 Subroutine MetricDefine3D(Dim,NC,NF,NP,IDS,FaceType,X,Y,Z,NFace_Cell,IFace_Cell,NConectCell,IConectCell,NConnectedPoints,IConnectedPoints,CellMetric,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NC,NF,NP
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:Dim)::FaceType,NFace_Cell,NConectCell,NConnectedPoints
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:100,1:Dim)::IConectCell,IConnectedPoints
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:3,1:3,1:Dim)::CellMetric,NodeMetric
!*********************************************************************************************
!Part 1:
 Call CellMetricDefination(Dim,NC,NF,NP,IDS,FaceType,X,Y,Z,NFace_Cell,IFace_Cell,CellMetric)
!Part 2:
 Call NodeMetricInterpolate(Dim,NC,NP,NConectCell,IConectCell,CellMetric,NodeMetric)
!Part 3:
 Call NodeMetricCorrecting3D(Dim,NC,NP,NConectCell,NConnectedPoints,IConnectedPoints,NodeMetric)
 
!*********************************************************************************************
 End
!###########################################################################################