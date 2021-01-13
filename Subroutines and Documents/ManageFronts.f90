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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,Corn,Neib,newQuad,current_Level
Intent(Inout)::FrontEdges,Fronts,States

Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4

Integer::Dim,NC,Fronts,current_Level,Elm,I,J,index
Integer,Dimension(1:2)::E
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::isFrontEdge
!===========================================================================================
do I=1,4
    
    Elm = Neib(NC,I)
    
    if(Elm /= 0) then
!Part 1:        
        if(Corn(Elm,4) == 0) then !------------------- if Elm is a triangle ----------------
            
            do J=1,3
                if(Neib(Elm,J) == NC) then
                    index = J
                    exit
                endif
            end do
            
            Fronts = Fronts + 1
            
            if(index == 1) then
                
                FrontEdges(Fronts,LeftVertex)=Corn(Elm,2)  ! Saving the LEFT corners of the discovered front edge
                FrontEdges(Fronts,RightVertex)=Corn(Elm,3) ! Saving the RIGHT corners of the discovered front edge 
                
            elseif(index == 2) then
            
                FrontEdges(Fronts,LeftVertex)=Corn(Elm,3)  ! Saving the LEFT corners of the discovered front edge
                FrontEdges(Fronts,RightVertex)=Corn(Elm,1) ! Saving the RIGHT corners of the discovered front edge 
                
            elseif(index == 3) then
            
                FrontEdges(Fronts,LeftVertex)=Corn(Elm,1)  ! Saving the LEFT corners of the discovered front edge
                FrontEdges(Fronts,RightVertex)=Corn(Elm,2) ! Saving the RIGHT corners of the discovered front edge
                
            endif
            
            FrontEdges(Fronts,Element) = Elm             ! Saving cell number of triangle elements having a Front Edge
            FrontEdges(Fronts,Level) = current_Level + 1 ! Saving Level of current Front Edge    
            
            print *,'Front Added!!! ',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)
            
        else !--------------------------- if Elm is a Quadrilateral -------------------------
!Part 2:        
            if(I /= 4) then
                E(1) = Corn(NC,I)
                E(2) = Corn(NC,I+1)
            else
                E(1) = Corn(NC,I)
                E(2) = Corn(NC,1)    
            endif
            
            do J=1,Fronts
		        if(States(J) /= Processed) then
			        if((FrontEdges(J,LeftVertex) == E(1) .And. FrontEdges(J,RightVertex) == E(2)) .Or. (FrontEdges(J,LeftVertex) == E(2) .And. FrontEdges(J,RightVertex) == E(1))) then
				        print *,'Front Processed: ',J,' : ',FrontEdges(J,LeftVertex),FrontEdges(J,RightVertex)
				        States(J) = Processed
			        endif
		        endif
	        end do
            
        endif
        
    else
!Part 3:        
        if(I /= 4) then
            E(1) = Corn(NC,I)
            E(2) = Corn(NC,I+1)
        else
            E(1) = Corn(NC,I)
            E(2) = Corn(NC,1)    
        endif
            
        do J=1,Fronts
		    if(States(J) /= Processed) then
			    if((FrontEdges(J,LeftVertex) == E(1) .And. FrontEdges(J,RightVertex) == E(2)) .Or. (FrontEdges(J,LeftVertex) == E(2) .And. FrontEdges(J,RightVertex) == E(1))) then
				    print *,'Front Processed: ',J,' : ',FrontEdges(J,LeftVertex),FrontEdges(J,RightVertex)
				    States(J) = Processed
			    endif
		    endif
	    end do
        
    endif
    
end do
!===========================================================================================
End Subroutine ManageFronts
!*********************************************************************************************
