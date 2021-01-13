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
 Subroutine Cons_New_Tri_Lawson(Dim,ME,Inode,NC,Corn,Neib,Nstack,Istack)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,ME,Inode
 Intent(Inout)::NC,Corn,Neib,Nstack,Istack

 Integer::Dim,NC,Inode,J,P1,P2,P3,N1,N2,N3,ME,Nstack
 Integer,Dimension(1:Dim,1:4)::Corn,Neib
 Integer,Dimension(1:Dim,1:2)::Istack
!********************************************************************************************* 
!part 1:
 P1=Corn(ME,1)
 P2=Corn(ME,2)
 P3=Corn(ME,3)

!part 2:
 N1=Neib (ME,1)
 N2=Neib (ME,2)
 N3=Neib (ME,3)

!part 3: 
 Corn(ME,1)=P1
 Corn(ME,2)=P2
 Corn(ME,3)=Inode

!part 4: 
 Corn(NC+1,1)=P2
 Corn(NC+1,2)=P3
 Corn(NC+1,3)=Inode

!part 5: 
 Corn(NC+2,1)=P3
 Corn(NC+2,2)=P1
 Corn(NC+2,3)=Inode

!part 6: 
 Neib (ME,1)=NC+1
 Neib (ME,2)=NC+2
 Neib (ME,3)=N3

!part 7: 
 Neib (NC+1,1)=NC+2
 Neib (NC+1,2)=ME
 Neib (NC+1,3)=N1

!part 8: 
 Neib (NC+2,1)=ME
 Neib (NC+2,2)=NC+1
 Neib (NC+2,3)=N2

!part 9:
 Do J=1,3
    If(N1/=0)Then
     IF(Neib(N1,J)==ME) Neib(N1,J)=NC+1
    Endif
	If(N2/=0)Then
     IF(Neib(N2,J)==ME) Neib(N2,J)=NC+2
    Endif
 End Do

!part 10: 
 Nstack=0
 If(N1/=0)Then
 Nstack=Nstack+1
 Istack(Nstack,1)=NC+1
 Istack(Nstack,2)=N1
 Endif

 If(N2/=0)Then
 Nstack=Nstack+1
 Istack(Nstack,1)=NC+2
 Istack(Nstack,2)=N2
 Endif

 If(N3/=0)Then
 Nstack=Nstack+1
 Istack(Nstack,1)=ME
 Istack(Nstack,2)=N3
 Endif

!part 11:
 NC=NC+2
!*********************************************************************************************
 End 
!###########################################################################################

