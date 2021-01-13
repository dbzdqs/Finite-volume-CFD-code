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
!// Date: Oct., 05, 2016                                                                   //!
!// Developed by: *//*-+/                       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Swap(Dim,ME,NE,Corn,Neib)
Implicit None
!===========================================================================================
 Intent(In   )::Dim,ME,NE
 Intent(Inout)::Corn,Neib
 
 Integer::Dim,J,I1,I2,J1,J2,ME,NE,P1,P2,P3,P4,N1,N2,N3,N4,Iface1,Iface2,Pi1,Pi2,Pi3,Pi4
 Integer,Dimension(1:Dim,1:4)::Corn,Neib
!===========================================================================================
!Part 1:
 Do J=1,3
	If( Neib(ME,J)==NE ) Iface1=J
	If( Neib(NE,J)==ME ) Iface2=J 
 End Do

!Part 2:
 If(    Iface1==1)Then
  I1=2
  I2=3
 Elseif(Iface1==2)Then
  I1=3
  I2=1
 Elseif(Iface1==3)Then
  I1=1
  I2=2
 Endif

!Part 3:
 P1=Corn(ME,I1)
 P2=Corn(ME,I2) 
 P3=Corn(ME,Iface1)
 P4=Corn(NE,Iface2)

!Part 4:
 Do J=1,3
	If( P1==Corn(NE,J) ) J1=J
	If( P2==Corn(NE,J) ) J2=J
 End Do

!Part 5:
 N1=Neib(ME,I1)
 N2=Neib(ME,I2)
 N3=Neib(NE,J1)
 N4=Neib(NE,J2)

!Part 6:                     
 Corn(ME,1)=P3
 Corn(ME,2)=P1
 Corn(ME,3)=P4
  
 Corn(NE,1)=P3
 Corn(NE,2)=P4
 Corn(NE,3)=P2

!Part 7:
 Neib (ME,1)=N4
 Neib (ME,2)=NE
 Neib (ME,3)=N2

 Neib (NE,1)=N3
 Neib (NE,2)=N1
 Neib (NE,3)=ME

!Part 8:
!Do J=1,4
!    If( N4/=0 .And. Neib(N4,J)==NE) Neib(N4,J)=ME
!    If( N1/=0 .And. Neib(N1,J)==ME) Neib(N1,J)=NE
!End Do
 
 !-------------------------------- Modifying N1 Neibours ----------------------
 if(N1/=0) then

	 if(Corn(N1,4)/=0) then

		do J=1,4
			if(Corn(N1,J)==P2) Pi2 = J
			if(Corn(N1,J)==P3) Pi3 = J
		end do

		if(Pi3 > Pi2) then
			if(Pi3 - Pi2 == 1) then
				Neib(N1,Pi2) = NE
			else
				Neib(N1,Pi3) = NE
			endif
		else
			if(Pi2 - Pi3 == 1) then
				Neib(N1,Pi3) = NE
			else
				Neib(N1,Pi2) = NE
			endif
		endif

	 else

		do J=1,3
			if(Corn(N1,J)/=P2 .And. Corn(N1,J)/=P3) then
				Neib(N1,J) = NE
				exit
			endif
		end do

	 endif
 
 endif

 !-------------------------------- Modifying N4 Neibours ----------------------
 if(N4/=0) then
 
	 if(Corn(N4,4)/=0) then
		
		do J=1,4
			if(Corn(N4,J)==P1) Pi1 = J
			if(Corn(N4,J)==P4) Pi4 = J
		end do

		if(Pi4 > Pi1) then
			if(Pi4 - Pi1 == 1) then
				Neib(N4,Pi1) = ME
			else
				Neib(N4,Pi4) = ME
			endif
		else
			if(Pi1 - Pi4 == 1) then
				Neib(N4,Pi4) = ME
			else
				Neib(N4,Pi1) = ME
			endif
		endif

	 else

		do J=1,3
			if(Corn(N4,J)/=P1 .And. Corn(N4,J)/=P4) then
				Neib(N4,J) = ME
				exit
			endif
		end do
 
	 endif 

endif
!===========================================================================================
 End Subroutine Swap
!*********************************************************************************************
