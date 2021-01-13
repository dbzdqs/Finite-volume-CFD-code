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
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine CheckThreeEdgeLoop(Dim,Corn,Neib,FrontEdges,States,Fronts,CF,isThreeEdgeLoop)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,FrontEdges,Fronts,CF
Intent(Out)::isThreeEdgeLoop
Intent(InOut)::States

Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Element=4

Integer::Dim,Fronts,CF,ME,LV,RV,NFE,NFI_L,NFI_R,NFEC,NF_Vertex,Count,I
Integer,Dimension(1:2)::Ei
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::isThreeEdgeLoop,isFrontEdge
!===========================================================================================
isThreeEdgeLoop = .True.
Count = 0

!Part 1:

LV = FrontEdges(CF,LeftVertex)
RV = FrontEdges(CF,RightVertex)
ME = FrontEdges(CF,Element)

!Part 2:

do I=1,3
    
    if(I /= 3) then
        Ei(1) = Corn(ME,I)
        Ei(2) = Corn(ME,I+1)
    else
        Ei(1) = Corn(ME,I)
        Ei(2) = Corn(ME,1)    
    endif
    
    if(isFrontEdge(Dim,Fronts,FrontEdges,States,Ei)) then
        Count = Count + 1    
    endif
    
end do

!Part 3:

if(Count == 3) then

    print *,'---------->>>>>>> Loop of Three Edges Eliminated! <<<<<<<----------'
    
    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,LV,RV,CF,1,Corn,Neib,States,NFE,NFI_L,NFEC,NF_Vertex)
    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,RV,LV,CF,2,Corn,Neib,States,NFE,NFI_R,NFEC,NF_Vertex)
    
    States(CF) = Processed
    States(NFI_L) = Processed
    States(NFI_R) = Processed
    
else
    isThreeEdgeLoop = .False.    
endif
   
!===========================================================================================
End Subroutine CheckThreeEdgeLoop
!*********************************************************************************************
