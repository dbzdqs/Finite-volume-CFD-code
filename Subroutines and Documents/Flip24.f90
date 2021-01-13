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
 Subroutine Flip24(Dim,NC,ME,NE,NP,Corn,Neib)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,ME,NE
 Intent(InOut)::NC,Corn,Neib
 
 Integer::Dim,I,J,I1,I2,ME,NE,P1,P2,P3,P4,N1,N2,N3,N4,Iface1,Iface2,NC,NP
 Integer,Dimension(1:Dim,1:4)::corn,Neib
!*********************************************************************************************
!Part 1:
 Do J=1,3
	If( Neib(ME,J)==NE ) Iface1=J
	If( Neib(NE,J)==ME ) Iface2=J 
 End Do

!Part 2:
 If(Iface1==1)Then
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
 P1=Corn(ME,Iface1)
 P2=Corn(NE,Iface2) 
 P3=Corn(ME,I1)
 P4=Corn(ME,I2)

!Part 4:
 Do J=1,3
	If( P3==Corn(ME,J) ) N1=Neib(ME,J)
	If( P4==Corn(ME,J) ) N2=Neib(ME,J)
	If( P4==Corn(NE,J) ) N3=Neib(NE,J)
	If( P3==Corn(NE,J) ) N4=Neib(NE,J)
 End Do

!Part 5:
 Corn(ME,1)=P2
 Corn(ME,2)=P4
 Corn(ME,3)=NP

 Neib(ME,1)=NC+1
 Neib(ME,2)=NE
 Neib(ME,3)=N4


 Corn(NE,1)=NP
 Corn(NE,2)=P3
 Corn(NE,3)=P2

 Neib(NE,1)=N3
 Neib(NE,2)=ME
 Neib(NE,3)=NC+2


 Corn(NC+1,1)=P4
 Corn(NC+1,2)=P1
 Corn(NC+1,3)=NP

 Neib(NC+1,1)=NC+2
 Neib(NC+1,2)=ME
 Neib(NC+1,3)=N1


 Corn(NC+2,1)=P1
 Corn(NC+2,2)=P3
 Corn(NC+2,3)=NP

 Neib(NC+2,1)=NE
 Neib(NC+2,2)=NC+1
 Neib(NC+2,3)=N2

!Part 6:
 Do J=1,3
    If( N2/=0 .And. Neib(N2,J)==ME ) Neib(N2,J)=NC+2
    If( N4/=0 .And. Neib(N4,J)==NE ) Neib(N4,J)=ME
    If( N1/=0 .And. Neib(N1,J)==ME ) Neib(N1,J)=NC+1
 End Do 

!Part 7:
 NC=NC+2
!*********************************************************************************************
 End 
!###########################################################################################