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
!// Developed by: M. Vakili, Computer Science, Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Read_STL(Dim,NP,NC,Corn,X,Y,Z)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim
 Intent (Out  )::NP,NC,Corn,X,Y,Z

 Integer::Dim,I,J,NP,NC,IP,IO,EOF
 Real*8::Rex,Rey,Rez,Eps
 Integer,Dimension(1:4,1:Dim)::Corn
 Real*8,Dimension(1:Dim,1:3)::Xd,Yd,Zd
 Real*8,Dimension(1:Dim)::X,Y,Z
 Character*200::TempStr,RowStr
!*********************************************************************************************
!Part 1:
 Open(1,File='STLIn.Stl')

!Part 2:
 Read(1,*)
 Read(1,*)

!Part 3:
 NC=0
 EOF=0

!Part 4:
 Do While(EOF/=1)

!Part 5:
    Read(1,*)
    NC=NC+1

!Part 6:
 	Do J=1,3
       Read(1,*) TempStr, Xd(NC,J), Yd(NC,J), Zd(NC,J)
    End Do

!Part 7:
    Read(1,*)
    Read(1,*)

!Part 8:
    Read(1,'(A200)',IOSTAT = IO) RowStr

!Part9   
    If(RowStr(1:8)=="endsolid") Then
     Read(1,'(A200)',IOSTAT = IO) RowStr
	 
	 If(RowStr(1:5)=="solid") Then
      Read(1,'(A200)',IOSTAT = IO) RowStr
     Else
       EOF=1
     End If
     
    End If

 End Do

!Part 10:
 Eps=0.00000001
 NP=0

!Part 11:
 Do I=1,NC
    Do J=1,3

	  !Part 12:
       Do IP=1,NP

          Rex=Dabs( Xd(I,J)-X(IP) )
          Rey=Dabs( Yd(I,J)-Y(IP) )
          Rez=Dabs( Zd(I,J)-Z(IP) )

          If( Rex<Eps .And. Rey<Eps .And. Rez<Eps )Then
           Corn(J,I)=IP
           Goto 10
          Endif

       End Do

	  !Part 13:
       NP=NP+1
       Corn(J,I)=NP

       X(Np)=Xd(I,J)
       Y(Np)=Yd(I,J)
       Z(Np)=Zd(I,J)

10  End Do
 End Do

 Close(1)
!*********************************************************************************************
 End Subroutine
!###########################################################################################
