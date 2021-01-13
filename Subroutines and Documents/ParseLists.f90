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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine ParseLists(Dim,NBC,BFP,States,Angles,current_Level,Corn,Neib,X,Y,Fronts,FrontEdges,NC,NP)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBC,BFP,current_Level
Intent(Inout)::Corn,Neib,X,Y,States,Angles,Fronts,FrontEdges,NC,NP

Integer,Parameter::Processed = -1
Integer,Parameter::STATUS_ZERO = 0
Integer,Parameter::STATUS_ONE = 1
Integer,Parameter::STATUS_TWO = 2
Integer,Parameter::STATUS_THREE = 3
Integer,Parameter::NO_CANDIDATE = 0
Integer,Parameter::LeftVertex = 1
Integer,Parameter::RightVertex = 2

Integer::Dim,CF,C,Fronts,NP,NC,NBC,current_Level,of_State,DF,I,J,K,L
Integer,Dimension(1:Dim)::States,DelayList,FList
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::FrontSelected,Done,Delayed
Real(8),Dimension(1:Dim)::X,Y
Real(8),Dimension(1:Dim,1:2)::Angles
!===========================================================================================
!Part 1:
FrontSelected = .False.
    
if(.Not. FrontSelected) then
	Call SelectCandidateFront(Dim,Corn,Neib,Fronts,FrontEdges,States,X,Y,current_Level,STATUS_THREE,CF,of_State)
	if(CF /= NO_CANDIDATE) then
		Call DefineSideEdges(Dim,NC,NP,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,Angles,X,Y,current_Level,CF,of_State)
        FrontSelected = .True.
    endif
endif

!Part 2:

if(.Not. FrontSelected) then
	Call SelectCandidateFront(Dim,Corn,Neib,Fronts,FrontEdges,States,X,Y,current_Level,STATUS_ONE,CF,of_State)
	if(CF /= NO_CANDIDATE) then
		Call DefineSideEdges(Dim,NC,NP,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,Angles,X,Y,current_Level,CF,of_State)
        FrontSelected = .True.
    endif
endif

!Part 3:

if(.Not. FrontSelected) then
	Call SelectCandidateFront(Dim,Corn,Neib,Fronts,FrontEdges,States,X,Y,current_Level,STATUS_TWO,CF,of_State)
	if(CF /= NO_CANDIDATE) then
		Call DefineSideEdges(Dim,NC,NP,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,Angles,X,Y,current_Level,CF,of_State)
        FrontSelected = .True.
    endif
endif

!Part 4:

if(.Not. FrontSelected) then
	Call SelectCandidateFront(Dim,Corn,Neib,Fronts,FrontEdges,States,X,Y,current_Level,STATUS_ZERO,CF,of_State)
	if(CF /= NO_CANDIDATE) then
		Call DefineSideEdges(Dim,NC,NP,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,Angles,X,Y,current_Level,CF,of_State)
        FrontSelected = .True.
    endif
endif
!===========================================================================================
End Subroutine ParseLists
!*********************************************************************************************
