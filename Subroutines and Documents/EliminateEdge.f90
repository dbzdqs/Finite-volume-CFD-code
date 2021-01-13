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
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
  Subroutine EliminateEdge(Dim,NContractedEdges,IContractedEdges,NEdgeOfCell,IEdgeOfCell,BeginOfReg,NR,NFR,NP,NF,NC,IDS,X,Y,NBP,IBP)
 Implicit None
!*********************************************************************************************
 Integer::Dim,I,J,F,S,K,NP,NF,NC,NR,P1,P2,E,EI,NContractedEdges,NConectCell,Heir,Dead,NegativeVol,Edgeexist,NBP,Point,Face
 Integer,Dimension(1:100)::NFR,BC,BeginOfReg
 Real(8),Dimension(1:Dim)::X,Y
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NEdgeOfCell
 Integer,Dimension(4,Dim)::IEdgeOfCell
 Integer,Dimension(1:Dim)::NConectEdge
 Integer,Dimension(1:Dim,1:100)::IConectEdge
 Integer,Dimension(1:50)::IConectCell
 Integer,Dimension(1:2,1:NF)::IContractedEdges
 Integer,Dimension(1:Dim)::IBP
 Logical::IsInvalidEdge=.FALSE.,repeated
 Integer,Dimension(1:4,1:Dim)::IPoint_Cell
 Integer,Dimension(1:Dim)::NPoint_Cell
 Real(8)::Area,DArea
 Integer::JJ,Cell
!**********************************************************************************************    
!Part 1:
 Do I=1,NContractedEdges      
     
   !Part 2:
    Heir = IContractedEdges(1,I)
    Dead = IContractedEdges(2,I)
        
   !Part 3:
    Call ConectedEdgeOfPointV2(Dim,NF,NP,BeginOfReg,NR,NFR,IDS,NConectEdge,IConectEdge)
     
   !Part 4:
    Edgeexist = -1
    Do J=1,NConectEdge(Dead)
       E    = IConectEdge(Dead,J)
       P1   = IDS(3,E)
       P2   = IDS(4,E)    
       IF((P1==Dead .and. P2==Heir) .or. (P2==Dead .and. P1==Heir))Then
        Edgeexist = 1
        EI = E
        exit
       End IF
    End Do
    If(Edgeexist==-1) Cycle
       
   !Part 5:
    Call ConectedCellOfPoint(Dim,Dead,IDS,NConectEdge,IConectEdge,NConectCell,IConectCell)
     
   !Part 6:
    Call VolumeCheck(Dim,Dead,heir,NConectCell,IConectCell,IEdgeOfCell,NEdgeOfCell,IDS,X,Y,NegativeVol)
       
   !Part 7:
    IF( NegativeVol==-1 )Then
     Call EdgeContraction2D(Dim,EI,Dead,Heir,NConectEdge,IConectEdge,NEdgeOfCell,IEdgeOfCell,NF,NC,IDS,NR,NFR,BeginOfReg)
     Call EdgeOfCellV3(Dim,NF,NC,BeginOfReg,NR,NFR,IDS,NEdgeOfCell,IEdgeOfCell) 
    Endif

 End Do
!*********************************************************************************************
 End 
!###########################################################################################