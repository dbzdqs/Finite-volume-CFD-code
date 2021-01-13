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
 Subroutine Equalize_Orient(Dim,FirstEO,NC,Corn,Neib)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NC,FirstEO
 Intent (Inout)::Corn,Neib

 Integer::Dim,I,J,J1,N1,N2,M1,M2,ME,NE,NC,Nstack,Dumy,FirstEO
 Integer,Dimension(1:4,1:Dim)::Neib,Corn
 Integer,Dimension(1:Dim)::Tri_Dirc
 Integer,Dimension(1:Dim,1:2)::Istack
!*********************************************************************************************
!Part 1:
 Do J=1,NC
    Tri_Dirc(J)=0
 End Do
 Tri_Dirc(FirstEO)=1

!Part 2:
 Nstack=0
 Do J1=1,3
    if(Neib(J1,FirstEO)==0)cycle
    Nstack=Nstack+1
    Istack(Nstack,1)=FirstEO
	Istack(Nstack,2)=Neib(J1,FirstEO)
 End Do

!Part 3:
 Do While(Nstack/=0)

   !Part 4:
    ME=Istack(Nstack,1)
	NE=Istack(Nstack,2)
    Nstack=Nstack-1

   !Part 5:
    if(Tri_Dirc(NE)==1)cycle

   !Part 6:
    Do J1=1,3
	   If(Neib(J1,ME)==NE)Then

	    If(J1==1)Then
	     M1=2
		 M2=3
	    Elseif(J1==2)Then
	     M1=3
		 M2=1
	    Elseif(J1==3)Then
	     M1=1
		 M2=2
		Endif

		Exit

	   Endif
    End Do

   !Part 7:
    Do J1=1,3
	   If(Neib(J1,NE)==ME)Then

	    If(J1==1)Then
	     N1=2
		 N2=3
	    Elseif(J1==2)Then
	     N1=3
		 N2=1
	    Elseif(J1==3)Then
	     N1=1
		 N2=2
		Endif

		Exit

	   Endif
    End Do

   !Part 8:
    M1=Corn(M1,ME)
	M2=Corn(M2,ME)
    N1=Corn(N1,NE)
	N2=Corn(N2,NE)

   !Part 9:
	If( N1==M1 .And. N2==M2 )Then
	 Dumy       = Corn(2,NE)
     Corn(2,NE) = Corn(3,NE)
	 Corn(3,NE) = Dumy

	 Dumy       = Neib(2,NE)
     Neib(2,NE) = Neib(3,NE)
	 Neib(3,NE) = Dumy
	Endif

   !Part 10:
	Tri_Dirc(NE)=1

   !Part 11:
    Do J1=1,3
	   Dumy=Neib(J1,NE)
	   If( Dumy/=0 .and. Tri_Dirc(Dumy)==0 )Then
        Nstack=Nstack+1
        Istack(Nstack,1)=NE
	    Istack(Nstack,2)=Dumy
	   Endif
    End Do

 End Do
!*********************************************************************************************
 End 
!###########################################################################################
