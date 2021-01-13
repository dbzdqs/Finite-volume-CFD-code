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
 Subroutine BE_Recovery1_2D(Dim,NP,NC,NBoundCrv,NFacCrv,BFacPt,Corn,Neib,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NBoundCrv
 Intent(Inout)::NP,NC,NFacCrv,BFacPt,Corn,Neib,X,Y

 Integer::Dim,I,J,J1,J2,J3,K,K1,K2,K3,I1,I2,P1,P2,P3,P4,ME,NE,SBE,NBoundCrv,NP,NC,Intersect
 Real(8)::X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros
 Integer,Dimension(1:100)::NFacCrv
 Integer,Dimension(1:Dim,1:2)::BFacPt
 Integer,Dimension(1:Dim,1:4)::Corn,Neib
 Real(8),Dimension(1:Dim)::X,Y
!*********************************************************************************************
!Part 1:
 SBE=0
 Do J=1,NBoundCrv
	SBE = SBE + NFacCrv(J)
 End Do

!Part 2:
 I=0
 Do J=1,NBoundCrv
       !pause
    Do While( I<NFacCrv(J) )
	   I=I+1  
       
      !Part 3:
       P1=BFacPt(I,1)
   	   P2=BFacPt(I,2)
 !print*,i,p1,p2
 !pause
      !Part 4:
 	   Do K1=1,NC
	      Do K2=1,3
	         If( Corn(K1,K2)==P1 )Then
		      Do K3=1,3
	             If( Corn(K1,K3)==P2 ) Goto 10
		      End Do
		     Endif
	      End Do
	   End Do
  print*,i
      !Part 5:
       X1=X(P1)
	   Y1=Y(P1)
	   X2=X(P2)
	   Y2=Y(P2)

      !Part 6:
	   Do K1=1,NC
	      Do K2=1,3
	         If( Corn(K1,K2)==P1 )Then

             !Part 7:	
 	          NE=Neib(K1,K2)
              ME=K1

             !Part 8:
	          If(K2==1)Then
	           P3=Corn(ME,2)
		       P4=Corn(ME,3)
	          Elseif(K2==2)Then
	           P3=Corn(ME,3)
	           P4=Corn(ME,1)
	          Elseif(K2==3)Then
	           P3=Corn(ME,1)
	           P4=Corn(ME,2)
	          Endif 

             !Part 9:
	          X3=X(P3)
	          Y3=Y(P3)
	          X4=X(P4)
	          Y4=Y(P4)
       
	         !Part 10: 
	          Call Cross_Lines(X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros,Intersect)
       
              if(Intersect==1) Goto 20

		     Endif
	      End Do
	   End Do

      !Part 11:
20	   NP=NP+1
	   X(NP)=Xcros
	   Y(NP)=Ycros

      !Part 12:
	   Call Flip24(Dim,NC,ME,NE,NP,Corn,Neib)

      !Part 13:
	   Do K=SBE,I+1,-1
	   	  BFacPt(K+1,1) = BFacPt(K,1)
	   	  BFacPt(K+1,2) = BFacPt(K,2)
	   End Do
 	   SBE=SBE+1
	   NFacCrv(J)=NFacCrv(J)+1

 	   BFacPt(I,1) = P1
 	   BFacPt(I,2) = NP

	   BFacPt(I+1,1) = NP
 	   BFacPt(I+1,2) = P2

10  End Do
 End Do
!*********************************************************************************************
 End 
!###########################################################################################
