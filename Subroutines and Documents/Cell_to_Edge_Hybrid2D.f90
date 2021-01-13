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
 Subroutine Cell_to_Edge_Hybrid2D(Dim,NC,Corn,Neib,CellType,NF,IDS)
 Implicit None
!*********************************************************************************************
 Intent (In   )                 ::  Dim,NC,Corn,Neib,CellType
 Intent (Out  )                 ::  NF,IDS

 Integer                        ::  Dim,i,I1,I2,J,J1,J2,ME,NE,P1,P2,NC,NF
 Integer,Dimension(1:4,1:Dim)   ::  Corn,Neib,Tneib
 Integer,Dimension(1:4,1:Dim)   ::  IDS
 Integer,Dimension(1:Dim)       ::  CellType
!*********************************************************************************************
!Part 1:
 Do i=1,NC
    Tneib(:,i) = Neib(:,i)
 End Do
 
!Part 2:
 NF=0

!Part 3:
 Do J=1,NC
   Do J1=1,CellType(J)
   
     !Part 4:
      If( Tneib(J1,J)==-1 ) Cycle

     !Part 5:
      If(J1==1)Then
	   P1=Corn(2,J)	 
	   P2=Corn(3,J) 
      Elseif(J1==2)Then
          If (CellType(J)==3) Then
            P1=Corn(3,J)
	        P2=Corn(1,J)
          Elseif(CellType(J)==4) Then
            P1=Corn(3,J)
            P2=Corn(4,J)
          End If
      Elseif(J1==3)Then
          If (CellType(J)==3) Then
            P1=Corn(1,J)
	        P2=Corn(2,J)
          Elseif(CellType(J)==4) Then
            P1=Corn(4,J)
            P2=Corn(1,J)
          End If
      Elseif(J1==4)Then
          P1=Corn(1,J)
          P2=Corn(2,J) 
      Endif

     !Part 6:
      ME=J
      NE=Tneib(J1,J)	


     !Part 7:
	  NF=NF+1
      Write(*,*) NF

	  IDS(1,NF)=ME
	  IDS(2,NF)=NE
	  IDS(3,NF)=P1
	  IDS(4,NF)=P2
	   
     !Part 8:
	  Do J2=1,CellType(J)
	     If( NE/=0 ) Then
          IF( Tneib(J2,NE)==ME ) Tneib(J2,NE)=-1
         End IF
	  End Do
	   
    End Do
 End Do 
!*********************************************************************************************
 End Subroutine Cell_to_Edge_Hybrid2D
!###########################################################################################