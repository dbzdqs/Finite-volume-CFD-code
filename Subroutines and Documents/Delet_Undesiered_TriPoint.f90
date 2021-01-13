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
 Subroutine Delet_Undesiered_TriPoint(Dim,NC,NP,NBC,NFacCrv,BFacPt,Corn,Neib,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NBC,NFacCrv,BFacPt
 Intent(Inout)::NC,NP,Corn,Neib,X,Y

 Integer::Dim,I,J,J1,JJ,N,NC,Temp,NBE,Desire,Intersect,P1,P2,P3,M,N1,N2,NP,NBC
 Real(8)::X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros
 Integer,Dimension(1:Dim,1:4)::Corn,Neib,Tempn
 Integer,Dimension(1:Dim,1:2)::BFacPt
 Integer,Dimension(1:100)::NFacCrv
 Real(8),Dimension(1:Dim)::X,Y
!*********************************************************************************************
!Part 1:
 NBE=0
 Do J=1,NBC
	NBE = NBE + NFacCrv(J)
 End Do

!Part 2:
 X3 = 2.0*X(NP)
 Y3 = 2.0*Y(NP)

!Part 3:
 N=0

!Part 4:
 Do I=1,NC
   ! Print*,I,NC

   !Part 5:
    P1=Corn(I,1)
	P2=Corn(I,2)
	P3=Corn(I,3)

   !Part 6:
	If( P1>NP-3 .Or. P2>NP-3 .Or. P3>NP-3 ) Cycle

   !Part 7:
    Desire=-1

   !Part 8:
    X4 = ( X(P1)+X(P2)+X(P3) )/3
    Y4 = ( Y(P1)+Y(P2)+Y(P3) )/3

   !Part 9:
    M =0
    Do J1=1,NBE

      !Part 10:
       N1=BFacPt(J1,1)
       N2=BFacPt(J1,2)

      !Part 11:
       X1=X(N1)
	   Y1=Y(N1)
	   X2=X(N2)
	   Y2=Y(N2)

      !Part 12:
       Call Cross_Lines(X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros,Intersect)
	   If(Intersect==1) M=M+1

    End Do

   !Part 13:
    If( (M/2)/=(M/2.0) ) Desire=1

   !Part 14:		                  
    If(Desire==1)Then

     N=N+1
     Do J=1,3
        Temp      = Corn(I,J)
        Corn(I,J) = Corn(N,J)
        Corn(N,J) = Temp

        Temp       = Neib(I,J)
        Neib(I,J)  = Neib(N,J)
        Neib(N,J)  = Temp
     End Do 

    !Part 15:
     Do J=1,NC
        Do JJ=1,3
		   Tempn(J,JJ)=0
           If( Neib(J,JJ)==N )Then
		    Neib(J,JJ)=I
		    Tempn(J,JJ)=1
		   End If
        End Do
     End Do
 
    !Part 16:
     Do J=1,NC
	    Do JJ=1,3
           If( Neib(J,JJ)==I .And. Tempn(J,JJ)/=1) Neib(J,JJ)=N
	    End Do
     End Do
		
	End If
 
 End Do

!Part 17:
 NC=N
 Do I=1,NC
   Do J=1,3
      If( Neib(I,J)>NC ) Neib(I,J)=0
   End Do
 End Do

!Part 18:
 NP=NP-3
!*********************************************************************************************
 End 
!###########################################################################################
