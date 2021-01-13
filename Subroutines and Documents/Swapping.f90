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
 Subroutine Swapping(Dim,ME,NE,Corn,Neib,Nstack,Istack)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,ME,NE
 Intent(Inout)::Corn,Neib,Nstack,Istack
 
 Integer::Dim,I,J,I1,I2,J1,J2,ME,NE,P1,P2,P3,P4,N1,N2,N3,N4,Nstack,Iface1,Iface2
 Integer,Dimension(1:Dim,1:4)::Corn,Neib
 Integer,Dimension(1:Dim,1:2)::Istack
!*********************************************************************************************
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
 Do J=1,3
    If( N4/=0 )Then
     IF(Neib(N4,J)==NE) Neib(N4,J)=ME
    endif
    If( N1/=0 )Then
     IF(Neib(N1,J)==ME) Neib(N1,J)=NE
    endif
 End Do 

!Part 9: 
 If(N1/=0)Then
  Nstack=Nstack+1
  Istack(Nstack,1)=NE
  Istack(Nstack,2)=N1
 End If
 
 If(N2/=0)Then
  Nstack=Nstack+1
  Istack(Nstack,1)=ME
  Istack(Nstack,2)=N2
 End If
 
 If(N3/=0)Then
  Nstack=Nstack+1
  Istack(Nstack,1)=NE
  Istack(Nstack,2)=N3
 End If

 If(N4/=0)Then
  Nstack=Nstack+1
  Istack(Nstack,1)=ME
  Istack(Nstack,2)=N4
 End If
!*********************************************************************************************
 End 
!###########################################################################################