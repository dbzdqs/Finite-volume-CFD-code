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
 Subroutine Find_Set_Elem_PointV3(Dim,NP,NC,Corn,N_Eset,I_Eset,N_Pset,I_Pset)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NC,Corn
 Intent(Out  )::N_Eset,I_Eset,N_Pset,I_Pset

 Integer::Dim,I,J,JJ,II,J1,J2,J3,N,P,P1,P2,ME,NP,NC
 Integer,Dimension(1:4,1:Dim)::Corn
 Integer,Dimension(1:Dim,1:50)::I_Pset,I_Eset
 Integer,Dimension(1:Dim)::N_Eset,N_Pset,NCorn
!*********************************************************************************************
!Part 1:
 Do J=1,NP    
    N_Eset(J) = 0 
 End Do

!Part 2:
 Do J=1,NC
    NCorn(J)=3
    IF( Corn(4,J)/=0 ) NCorn(J)=4
 End Do

!Part 3:
 Do J=1,NC

   !Part 4:
    Do J1=1,NCorn(J)  
       P=Corn(J1,J)
        
      !Part 5: 
       N_Eset(P) = N_Eset(P) + 1
       I_Eset(P,N_Eset(P)) = J
       
    End Do

 End Do

!Part 6:
 Do J=1,NP

   !Part 7:
    N=0

   !Part 8:
    Do J1=1,N_Eset(J)
       ME=I_Eset(J,J1)

      !Part 9:
       Do J2=1,NCorn(ME)
	      P=Corn(J2,ME)
	
         !Part 10:
	      If( P==J )Cycle
	
         !Part 11:
          Do J3=1,N
             If( P==I_Pset(J,J3) ) Goto 10
		  End do
	
         !Part 12:
          N=N+1
          I_Pset(J,N) = P

10     End Do

    End Do
	
   !Part 13:
	N_Pset(J) = N

 End Do
!*********************************************************************************************
 End 
!###########################################################################################
