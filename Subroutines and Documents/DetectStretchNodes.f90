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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DetectStretchNodes(Dim,NP,NodeMetric,StretchNode,ElipsoidRatio)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,EDir
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Logical,Dimension(1:Dim)::StretchNode
 Real(8)::MaxVal,MinVal,AspectRatio
 Real(8),Dimension(1:Dim)::ElipsoidRatio
!*********************************************************************************************
!Part 1:
 ElipsoidRatio(1:NP) = 0.0
 StretchNode(1:NP)   = .FALSE.
 Do I=1,NP
   !Part 2:
    Call CalcMinMaxLocalSize(NodeMetric(:,:,I),MinVal,MaxVal,EDir)
    AspectRatio = MaxVal/MinVal
    
   !Part 3:
    If(AspectRatio >= 1.30) StretchNode(I) = .TRUE.
    ElipsoidRatio(I) = AspectRatio
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################