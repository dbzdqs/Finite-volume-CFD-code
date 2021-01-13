!DDDDDDDDDDDDDDDDDDDDDDDDDDDAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNSSSSSSSSSSSSSSSSSSSSSSSSSS
!//             /////////////       ////////////////    ////////     //////    ////////////////        //!
!//             /////////////      ////////////////    //////////   //////    ////////////////         //!
!//            /////    /////     //////    //////    //////////// //////    /////                     //!
!//           /////    //////    ////////////////    ///////////////////    ////////////////           //!
!//          /////    //////    ////////////////    ////// ////////////               /////            //!
!//         ///////////////    //////    //////    //////   //////////    ////////////////             //!
!//       ///////////////     //////    //////    //////     ////////    ////////////////              //!
!//          Developer            Assistant    in      Numerical             Sciences                  //!
!//----------------------------------------------------------------------------------------------------//!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace Department, Amirkabir University of Technology           //!
!// Chief Developer: N. msnkre, Aerospace eng., Amirkabir University of Technology                     //!
!// Date: October, 14, 2013                                                                            //!
!//                                                                                                    //!
!// The Program is Available Through the Website: www.DANS.ir                                          //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                               //!
!//----------------------------------------------------------------------------------------------------//!
!// Description:                                                                                       //!
!// The main objective of this code is the generation of computational mesh by Q-Morph method. In the  //!
!// first step of this algorithm is converting triangular mesh to quadrilateral mesh. In the next      //!
!// step, quality of mesh increased by using blended algorithms that improve topology of the mesh and  //!
!// finally  by using blended algorithms which is based on systematic changing the location of the     //!
!// nodes, quality of elements improved .                                                              //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
PROGRAM QMorph
Implicit None
!===============================
Integer,Parameter::Dim=80000
!===============================
Integer::I,J,Step
Integer::Elm
Integer::TEC
Integer::QEC
Integer,Dimension(1:1000)::TElms
Integer,Dimension(1:1000)::QElms
Logical::ElementInversionOccured
Logical::QuadIsInverted
Integer::NP  !Number of Existing Points
Integer::NC !Number of Cells of mesh
Integer::NBC !Number of Boundary Curves
Integer::NBE
Integer::Fronts
Integer::current_Level
Integer::Num_Of_Loops
Integer,Dimension(1:10)::NFC
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::FrontEdges
Integer,Dimension(1:Dim,1:4)::Corn !Corners point index of Each Cell 
Integer,Dimension(1:Dim,1:4)::Neib !Neighboring Cell Index of Each Elements 
Logical::AllFrontsProcessed
Logical::Continuity
Logical::CurrentLevelCompleted
Logical::NumOfActiveFrontsIsOdd
Logical::done
Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
Real(8),Dimension(1:Dim,1:2)::Angles
!*********************************************************************************************************
!Part 1:
I=0
NBE=0
Fronts=0
current_Level=1

Call Initialize(Dim,FrontEdges,States)
Call Read_2DMeshC_TriToQuad(Dim,NP,NC,NBC,NFC,BFP,Corn,Neib,X,Y)
Call Get_FrontEdges(Dim,Corn,Neib,NC,FrontEdges,Fronts,current_Level)
Call CheckContinuity(Dim,Fronts,FrontEdges,States,Num_Of_Loops,Continuity)

if(Continuity) then

!Part 2:

	Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

	do J=1,NBC
		NBE = NBE + NFC(J)
	end do

!Part 3:

	do
		if(CurrentLevelCompleted(Dim,FrontEdges,Fronts,States,current_Level)) then

			print*,'*****>>>Level:',current_Level,' Completed!!! <<<*****'

			current_Level = current_Level + 1

		endif

		if(AllFrontsProcessed(Dim,Fronts,States)) then
			exit
        endif

		print *,'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ROUND:',I

		Call ParseLists(Dim,NBE,BFP,States,Angles,current_Level,Corn,Neib,X,Y,Fronts,FrontEdges,NC,NP)

		I = I + 1

	end do
	
else
	print *,'>>>>>>>>>> Unfortunately Given Mesh Is NOT Standard <<<<<<<<<<'
endif

!Part 4:

if(AllFrontsProcessed(Dim,Fronts,States)) then
        
	Call TopologicalCleanUp(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y)
	Call GlobalSmooth(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP)
    
    do step=1,50
        
        I = 0
        print*,'Step: ',step
       
        Call AddElements(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y)

        Call TopologicalCleanUp(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y)       
        Call GlobalSmooth(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP)
 
        do J=1,NC
            if(Corn(J,1)/=Corn(J,2) .And. Corn(J,4) /= 0) then
                Call SizeCleanup(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y,J,done)
                if(done) then
                    I = I + 1
                    print*,J,':SizeCleanup'    
                endif
            endif
        end do

        Call TopologicalCleanUp(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y)        
        Call GlobalSmooth(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP)
         
        if(I == 0) exit
        
    end do
    
endif

Call CheckContinuity(Dim,Fronts,FrontEdges,States,Num_Of_Loops,Continuity)

print*,'Current Level: ',current_Level

print*,'Fronts: ',Fronts

!Part 5:
 Call Write_2DMeshC_TriToQuad(Dim,NP,NC,NBC,NFC,BFP,Corn,Neib,X,Y)
!********************************************************************************************************
 END
!########################################################################################################