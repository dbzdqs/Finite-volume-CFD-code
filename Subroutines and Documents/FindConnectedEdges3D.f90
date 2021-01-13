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
 Subroutine FindConnectedEdges3D(Dim,NF,NP,IDS,FaceType,NConnectedPoints,IConnectedPoints)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,I,J,k,P1,P2,Face,Q,P
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints
 Integer,Dimension(1:Dim)::NConnectedPoints
 Integer,Dimension(1:4)::FacePoints
 Logical::exists
!*********************************************************************************************
!Part 1:
 NConnectedPoints(:) = 0
 Do I=1,NF
     
 !Part 2:     
    Face = I
    FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
       
   !Part 3: 
    Do J=1,FaceType(Face)
       P = FacePoints(J)
       Q = FacePoints(Mod(J,FaceType(Face)) + 1)
           
      !Part 4: 
       exists = .FALSE.
       Do K=1,NConnectedPoints(P)
          If(IConnectedPoints(K,P)==Q)Then
           exists = .TRUE.
           exit
          Endif
       EndDo
       If(.NOT. exists)Then
        NConnectedPoints(P) = NConnectedPoints(P)+1
        IConnectedPoints(NConnectedPoints(P),P) = Q
       Endif
       
       exists = .FALSE.
       Do K=1,NConnectedPoints(Q)
          If(IConnectedPoints(K,Q)==P)Then
           exists = .TRUE.
           exit
          Endif
       EndDo
       If(.NOT. exists)Then
        NConnectedPoints(Q) = NConnectedPoints(Q)+1
        IConnectedPoints(NConnectedPoints(Q),Q) = P
       Endif
           
    EndDo
 EndDo
!*********************************************************************************************
 End
!###########################################################################################