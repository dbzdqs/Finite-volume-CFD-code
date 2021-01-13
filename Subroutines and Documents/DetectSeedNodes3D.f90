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
 Subroutine DetectSeedNodes3D(Dim,NP,NF,BeginOfReg,NFR,IDS,FaceType,LayerIndex,NSeedNodes,ISeedNodes)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NSeedNodes,I,J,K,P1,P2,Face
 Integer,Dimension(1:Dim)::LayerIndex,ISeedNodes,FaceType
 Logical::existsP1,existsP2
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:100)::BeginOfReg,NFR
 Integer,Dimension(1:4)::FacePoints
!*********************************************************************************************
!Part 1:
 NSeedNodes = 0
 Do I=BeginOfReg(2),(BeginOfReg(2)+NFR(2)-1)
    Face = I
    FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
    
   !Part 2:
    Do J=1,FaceType(Face)    
       P1 = FacePoints(J)
       P2 = FacePoints(Mod(J,FaceType(Face)) + 1)
          
      !Part 3:
       existsP1 = .FALSE.
       existsP2 = .FALSE.
       Do K=1,NSeedNodes
          If(ISeedNodes(K)==P1) existsP1 = .TRUE.
          If(ISeedNodes(K)==P2) existsP2 = .TRUE.
       EndDo
          
      !Part 4:
       If(.NOT. existsP1)Then
        NSeedNodes = NSeedNodes + 1
        ISeedNodes(NSeedNodes) = P1
        LayerIndex(P1) = 1
       Endif
          
       If(.NOT. existsP2)Then
        NSeedNodes = NSeedNodes + 1
        ISeedNodes(NSeedNodes) = P2
        LayerIndex(P2) = 1
       Endif
          
    EndDo
       
 EndDo
!*********************************************************************************************
 End
!###########################################################################################