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
 Subroutine Smart_Lap_Smooth(Dim,NP,Corn,N_Eset,I_Eset,N_Pset,I_Pset,X,Y)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NP,Corn,N_Eset,I_Eset,N_Pset,I_Pset
 Intent (InOut)::X,Y
 
 Real(8)::SumX,SumY,Xhat,Yhat,Min,Qual
 Integer::Dim,I,J,JJ,NP,P1,P2,P3,P_Temp,NEset,NPset
 Integer,Dimension(1:Dim,1:4)::Corn
 Integer,Dimension(1:Dim,1:50)::I_Eset,I_Pset
 Integer,Dimension(1:Dim)::N_Eset,N_Pset
 Real(8),Dimension(1:Dim)::X,Y,Quality
 Real(8),Dimension(1:50)::Temp_Qual
!*********************************************************************************************
!Part 1:
 Do I=1,NP

   !Part 2:
	NEset=N_Eset(I)
	NPset=N_Pset(I)

   !Part 4:
    If(NEset==0)Cycle


   !Part 5:   
    Min=100000000.0
    Do J=1,NEset
       JJ=I_Eset(I,J)

	   P1 =Corn(JJ,1)
       P2 =Corn(JJ,2)
       P3 =Corn(JJ,3)
       Call Quality_AreaLen(Dim,P1,P2,P3,X,Y,Qual)
       
	   Quality(JJ) = Qual

	   If( Qual<Min ) Min = Qual
    End Do

   !Part 4:
	SumX=0.0
	SumY=0.0
    Do J=1,NPset
       JJ=I_Pset(I,J)
 	   SumX = SumX + X(JJ)
	   SumY = SumY + Y(JJ)
    End Do
 	Xhat = SumX / NPset
	Yhat = SumY / NPset

   !Part 6:
    Do J=1,NEset
       JJ = I_Eset(I,J)
       
      !Part 7:
       P1 = Corn(JJ,1)
       P2 = Corn(JJ,2)
       P3 = Corn(JJ,3)

      !Part 8:
       P_Temp=NP+1
       X(P_Temp) = Xhat
       Y(P_Temp) = Yhat

      !Part 9:
       If(P1==I)P1=P_Temp       
       If(P2==I)P2=P_Temp       
       If(P3==I)P3=P_Temp

      !Part 10:
       Call Quality_AreaLen(Dim,P1,P2,P3,X,Y,Qual)
	   Temp_Qual(J) = Qual

	   If( Qual<Min ) Goto 10
    End Do
   
   !Part 11:
    X(I) = Xhat
    Y(I) = Yhat
    
   !Part 12:
    Do J=1,NEset
       JJ=I_Eset(I,J)
	   Quality(JJ) = Temp_Qual(J)
    End Do

10 End Do
!*********************************************************************************************
 End
!###########################################################################################

 