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
 Subroutine CellMetricDefination(Dim,NC,NF,NP,IDS,FaceType,X,Y,Z,NFace_Cell,IFace_Cell,CellMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NC,NF,NP,I,J,K,S,P1,P2,NCalcedEdges,Face,tempInt
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell,FaceType
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:3,1:3,1:Dim)::CellMetric
 Integer,Dimension(1:12,1:2)::CalcedEdges
 Logical::exists
 Integer,Dimension(1:2)::MDim1,MDim2
 Real(8),Dimension(1:3,1)::TNodeCP1P2
 Real(8),Dimension(1,1:3)::NodeCP1P2
 Real(8),Dimension(1:3,1:3)::TempMatrix
 Integer,Dimension(1:4)::FacePoints
!*********************************************************************************************
!Part 1:
 Do I=1,NC
     
    CellMetric(:,:,I) = 0
    NCalcedEdges      = 0
    Do J=1,NFace_Cell(I)
       
       Face = IFace_Cell(J,I)
       FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
        
      !Part 2: 
       Do K=1,FaceType(Face)
           P1 = FacePoints(K)
           P2 = FacePoints(Mod(K,FaceType(Face)) + 1)
           
           If(IDS(2,Face)==I)Then
               tempInt=P1
               P1=P2
               P2=tempInt
           Endif 
           
          !Part 3:
           exists = .FALSE.
           Do S=1,NCalcedEdges
              If((CalcedEdges(S,1)==P1 .AND. CalcedEdges(S,2)==P2) .OR. (CalcedEdges(S,1)==P2 .AND. CalcedEdges(S,2)==P1))Then
               exists = .TRUE.
               exit
              Endif
           EndDo  
           If(exists) Cycle
       
          !Part 4:
           NodeCP1P2(:,:)  = 0
           TNodeCP1P2(:,:) = 0
           NodeCP1P2(1,1)  = X(P2)-X(P1)
           NodeCP1P2(1,2)  = Y(P2)-Y(P1)
           NodeCP1P2(1,3)  = Z(P2)-Z(P1)
           TNodeCP1P2(1,1) = NodeCP1P2(1,1)
           TNodeCP1P2(2,1) = NodeCP1P2(1,2)
           TNodeCP1P2(3,1) = NodeCP1P2(1,3)
      
           NCalcedEdges =  NCalcedEdges + 1
           CalcedEdges(NCalcedEdges,1) = P1
           CalcedEdges(NCalcedEdges,2) = P2
                 
           MDim1(1) = 3
           MDim1(2) = 1     
           MDim2(1) = 1
           MDim2(2) = 3     
           TempMatrix(:,:) = 0
           Call MulMatrix(NodeCP1P2,MDim1,TNodeCP1P2,MDim2,TempMatrix)
       
           CellMetric(:,:,I) = CellMetric(:,:,I) + TempMatrix(:,:)
        EndDo    
    Enddo
     
   !Part 5:
    TempMatrix(:,:) = 0
    Call CalcMatrixInverseV2(CellMetric(:,:,I),TempMatrix(:,:),3)
    CellMetric(:,:,I)     = TempMatrix(:,:)
    CellMetric(1:3,1:3,I) = ((2.0)*CellMetric(1:3,1:3,I))
      
 EndDo
!*********************************************************************************************
 End
!###########################################################################################