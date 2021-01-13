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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Find_Neib(Dim,NC,Corn,Neib)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,Corn
 Intent(Out  )::Neib

 Integer::Dim,I,II,II1,II2,J,JJ,JJ1,JJ2,N1,N2,M1,M2,NC
 Integer,Dimension(1:4,1:Dim)::Corn,Neib
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    Neib(1,i)=0
    Neib(2,i)=0
    Neib(3,i)=0
    Neib(4,i)=0
 End do

!Part 2:
 Do I=1,NC
    Do II=1,3

	  !Part 3:
	   If(Neib(II,I)/=0)Cycle
		   
	  !Part 4:
	   If(II==1)Then
	    II1=2
		II2=3
	   Elseif(II==2)Then
	    II1=3
	    II2=1
	   Elseif(II==3)Then
	    II1=1
	    II2=2
	   Endif 
	   
	  !Part 5:
	   M1=Corn(II1,I)
	   M2=Corn(II2,I)
	   
	  !Part 6:
       Do J=1,NC
	   
	     !Part 7:
	      If(I==J)Cycle
	   
	     !Part 8:
	      Do JJ=1,3
	         If(Neib(JJ,J)/=0)Cycle
	  	   
	        !Part 9:
	         If(JJ==1)Then
	          JJ1=2
		      JJ2=3
	         Elseif(JJ==2)Then
	          JJ1=3
	          JJ2=1
	         Elseif(JJ==3)Then
	          JJ1=1
	          JJ2=2
	         Endif 
	   
	        !Part 10:
	         N1=Corn(JJ1,J)
	         N2=Corn(JJ2,J)
	   
	        !Part 11:
             If( (N1==M1 .And. N2==M2) .Or. (N1==M2 .And. N2==M1) )Then
	          Neib(II,I)=J     	   
              Neib(JJ,J)=I
	          Goto 10
	         Endif

          End Do	

       End Do

10  End Do
 End Do

!*********************************************************************************************
 End 
!###########################################################################################
