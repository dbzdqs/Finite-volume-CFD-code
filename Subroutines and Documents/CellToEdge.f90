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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CellToEdge(Dim,NC,Corn,Neib,NF,IDS)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,Corn,Neib,NC
 Intent (Out  )::IDS,NF

 Integer::Dim,I1,I2,J,J1,J2,ME,NE,P1,P2,NC,NF
 Integer,Dimension(1:Dim,1:4)::Corn,Neib,Tneib
 Integer,Dimension(1:4,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 Tneib=Neib

!Part 2:
 NF=0

!Part 3:
 Do J=1,NC
   Do J1=1,3
   
     !Part 4:
      If( Tneib(J,J1)==-1 ) Cycle

     !Part 5:
      If(J1==1)Then
	   P1=Corn(J,2)	 
	   P2=Corn(J,3) 
      Elseif(J1==2)Then
	   P1=Corn(J,3)	 
	   P2=Corn(J,1)    
	  Elseif(J1==3)Then
	   P1=Corn(J,1)	 
	   P2=Corn(J,2)
      Endif

     !Part 6:
      ME=J
      NE=Tneib(J,J1)	


     !Part 7:
	  NF=NF+1

	  IDS(1,NF)=ME
	  IDS(2,NF)=NE
	  IDS(3,NF)=P1
	  IDS(4,NF)=P2
	   
     !Part 8:
	  Do J2=1,3
	     If( NE/=0 ) Then
          IF( Tneib(NE,J2)==ME ) Tneib(NE,J2)=-1
         End IF
	  End Do
	   
    End Do
 End Do 
!*********************************************************************************************
 End 
!###########################################################################################
