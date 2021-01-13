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
!// Date: Nov., 15, 2014                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine TopologicalCleanUp(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP
Intent(InOut)::NC,NP,Corn,Neib,X,Y

Integer::Dim,NC,NP,NBE,I,J
Integer,Dimension(1:Dim)::IE
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::CornerFound,ElementFound,CleanUpPerformed,done,ElementInversionOccured
Real(8)::x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
J = 0
do

	do I=1,NC

		if(Corn(I,4)/=0 .And. Corn(I,1)/=Corn(I,2)) then
!Part 1:
            
			Call NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,I,done)
			
			if(done) then
			
				print *,'I: ',I,' exit! Node'

				exit

			endif

!Part 2:
			Call ElementElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,I,done)
			
			if(done) then
			
				print *,'I: ',I,' exit! Element'
			
				exit
			
			endif
!Part 3:
			Call SegmentElimination(Dim,Corn,Neib,X,Y,NBE,BFP,I,done)
				
			if(done) then
				
                J = J + 1
                
				print *,'I: ',I,' SegmentElimination',J

				exit
			
			endif
!Part 4:
			Call ThreeEdgedNodesDividedByFourEdgedNode(Dim,Corn,Neib,X,Y,NBE,BFP,I,done)
					
			if(done) then
			
				print *,'I: ',I,' ThreeEdgedNodesDividedByFourEdgedNode'

				exit

			endif
!Part 5:
			Call ThreeFiveOppoThree(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,I,done)
				
			if(done) then
			
				print *,'I: ',I,' ThreeFiveOppoThree'
			
				exit
			
			endif
!Part 6:
			Call FiveThreeOppoFive(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,I,done)

			if(done) then
			
				print *,'I: ',I,' FiveThreeOppoFive'
			
				exit
			
			endif
!Part 7:
			Call TwoCollapsesOperation(Dim,Corn,Neib,X,Y,NBE,BFP,I,done)

			if(done) then
			    
                
				print *,'I: ',I,' TwoCollapsesOperation'
			
				exit
			
			endif
!Part 8:
			Call TwoFiveOppoThree(Dim,Corn,Neib,X,Y,NBE,BFP,I,done)

			if(done) then
			
				print *,'I: ',I,' TwoFiveOppoThree'
			
				exit
			
			endif
!Part 9:
			Call ThreeFiveOppoFiveThree(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,I,done)

			if(done) then
			    
                !J = J + 1
                
				print *,'I: ',I,' ThreeFiveOppoFiveThree'
			
				exit
			
			endif
!Part 10:
			Call ThreeFiveFive(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,I,done)

			if(done) then
			
				print *,'I: ',I,' ThreeFiveFive'
			
				exit
			
            endif
			
		endif	

	end do

	if(I > NC) then
        print*,'--- End of TopologicalCleanup ---'
        exit
    endif

end do

!===========================================================================================
End Subroutine TopologicalCleanUp
!*********************************************************************************************
